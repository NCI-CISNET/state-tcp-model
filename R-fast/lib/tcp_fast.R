# Fast, output-equivalent implementation of the cigarette-tax model in R/02_model_functions.R
# (runstates, generate_prevs, calculate_mort, generate_TCPoutput, tax_effectCalculation). Every
# arithmetic operation is performed in the same order as the original, so results are
# bit-identical; tools/check_equivalence.R checks this against the original code.
#
# How it differs from the original (none of these change any output):
#   - State inputs, mortality reshaping and the baseline (no-policy) run are computed once per
#     state and kept (get_state), instead of being reloaded and rerun for every scenario and policy
#     year. In the original, every runstates() call reran the baseline.
#   - The cohort model's loops over cohorts, ages and calendar years are replaced by vectorised
#     matrix operations (cohort_prevs, sad_ap), or by the compiled kernel in cohort_kernel.cpp
#     when load_kernel() succeeds. A runstates(full = FALSE) call takes ~6 ms against 5-8 s for
#     the original runstates() (measured by tools/check_equivalence.R).
#   - The cohort model no longer builds m_output (a character matrix nothing reads); it is NULL.
#   - runstates(full = FALSE) skips everything the TCP tool files do not use: the former-smoker
#     prevalence arrays, the never-smoker age-period matrix, df_mort.outputs and
#     df_CSprevs.by.state.
#   - tcp_tables() builds the same tables as generate_TCPoutput() for many scenarios of one state:
#     the row layout and baseline columns, which do not depend on the scenario, are computed once
#     per state and policy year.
#   - generate_TCPoutput() only creates output directories when write_files = TRUE.

startbc <- 1908
endbc <- 2100
endyear <- 2200
cohyears <- endbc - startbc + 1
v_calyears <- 1:cohyears
v_stdbirths <- rep(1000000, times = 193)

statepop <- usmap::statepop

fips_to_abbr <- function(fips2) {
  fips2 <- sprintf("%02d", as.integer(fips2))
  abbr <- unname(setNames(statepop$abbr, statepop$fips)[fips2])
  if (is.na(abbr)) stop("Unknown state FIPS: ", fips2)
  abbr
}

abbr_to_fips <- function(abbr) {
  fips <- unname(setNames(statepop$fips, statepop$abbr)[toupper(abbr)])
  if (anyNA(fips)) stop("Unknown state abbreviation: ", paste(abbr[is.na(fips)], collapse = ", "))
  fips
}

# Cell (age r, cohort b) of an age-cohort matrix sits at (r, b + r - 1) of the age-period matrix;
# .ac_lin and .ap_lin are those cells' linear positions in each.
.nc <- cohyears
.r <- rep(1:100, times = .nc)
.b <- rep(1:.nc, each = 100)
.ap_lin <- .r + (.b + .r - 2) * 100
.ac_lin <- .r + (.b - 1) * 100
.ap_lin3 <- rep(.ap_lin, 40) + rep((0:39) * 100 * (.nc + 100), each = length(.ap_lin))

# Index sets for ap_to_ac() on the 100 x 193 SAD/YLL matrices.
.ac_c <- rep(1:.nc, each = 100)
.ac_ok <- .ac_c + .r - 1 <= .nc
.ac_tail <- !.ac_ok & .r >= 3
.ac_ok_dst <- (.r + (.ac_c - 1) * 100)[.ac_ok]
.ac_ok_src <- (.r + (.ac_c + .r - 2) * 100)[.ac_ok]
.ac_tail_dst <- (.r + (.ac_c - 1) * 100)[.ac_tail]
.ac_tail_src <- (.r + (.r - 3) * 100)[.ac_tail]

# Former-smoker compartments are held as a (cohort x years-since-quit) by age matrix in the
# R recurrence so each age step reads and writes one contiguous column.
.fs_shift <- 1:(38 * .nc)
.fs_39 <- 38 * .nc + 1:.nc
.fs_40 <- 39 * .nc + 1:.nc


#===============================================================================
# Optional compiled kernel
#===============================================================================

.kernel <- new.env()

load_kernel <- function(cpp_path, cache_dir = tools::R_user_dir("tcp_fast", "cache")) {
  ok <- tryCatch({
    old <- Sys.getenv("PKG_CXXFLAGS")
    Sys.setenv(PKG_CXXFLAGS = paste(old, "-ffp-contract=off"))
    on.exit(Sys.setenv(PKG_CXXFLAGS = old), add = TRUE)
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    Rcpp::sourceCpp(cpp_path, env = .kernel, cacheDir = cache_dir)
    TRUE
  }, error = function(e) {
    message("cohort_kernel.cpp could not be compiled, using the R implementation: ", conditionMessage(e))
    FALSE
  })
  invisible(ok)
}

use_kernel <- function() {
  isTRUE(getOption("tcp_fast.kernel", TRUE)) && exists("cohort_kernel", envir = .kernel, inherits = FALSE)
}


#===============================================================================
# State inputs and baseline, cached per process
#===============================================================================

.cache <- new.env()

# Survival probabilities (1 - death probability) by age and cohort. surv*_t are laid out for the
# R recurrence in cohort_prevs(); survNS_ac and survCS_ac are age x cohort for the kernel, and
# survYSQ_jt is years since quit x age x cohort, the kernel's layout for former smokers.
prepare_mortality <- function(p_mortNS_AC, p_mortCS_AC, p_mortYSQ_AC) {
  cohs <- startbc:endbc
  mY <- array(0, c(100, .nc, 40))
  rY <- split(seq_len(dim(p_mortYSQ_AC)[1]), p_mortYSQ_AC[, 1, 1])
  rNS <- split(seq_len(nrow(p_mortNS_AC)), p_mortNS_AC[, 1])
  rCS <- split(seq_len(nrow(p_mortCS_AC)), p_mortCS_AC[, 1])
  for (k in seq_along(cohs)) mY[, k, ] <- p_mortYSQ_AC[rY[[as.character(cohs[k])]], 5, ]
  mNS <- vapply(cohs, function(k) p_mortNS_AC[rNS[[as.character(k)]], 5], numeric(100))
  mCS <- vapply(cohs, function(k) p_mortCS_AC[rCS[[as.character(k)]], 5], numeric(100))
  sY <- 1 - mY
  list(survNS_ac = 1 - mNS, survCS_ac = 1 - mCS, survYSQ_jt = aperm(sY, c(3, 1, 2)),
       survNS_t = t(1 - mNS), survCS_t = t(1 - mCS),
       survYSQ_t = matrix(aperm(sY, c(2, 3, 1)), .nc * 40, 100))
}

prepare_mort_ap <- function(m_p_mortNS_AP, m_p_mortCS_AP, a_p_mortYSQ_AP, df_census_data) {
  m_popdist <- as.matrix(df_census_data)
  colnames(m_popdist) <- as.character(startbc:endbc)
  dY <- array(0, dim(a_p_mortYSQ_AP))
  for (j in 1:40) dY[, , j] <- a_p_mortYSQ_AP[, , j] - m_p_mortNS_AP
  list(popdist = m_popdist, dCS = m_p_mortCS_AP - m_p_mortNS_AP, dYSQ = dY,
       dYSQ_jt = aperm(dY, c(3, 1, 2)))  # years since quit x age x year, for the kernel
}

get_state <- function(fipscode, data_dir) {
  key <- paste(data_dir, fipscode)
  if (identical(.cache$key, key)) return(.cache$state)

  e <- new.env()
  load(file.path(data_dir, "state_inputs", "mort_rates", paste0("p.mort_", fipscode, ".RData")), envir = e)
  load(file.path(data_dir, "state_inputs", paste0("smk_", fipscode, ".RData")), envir = e)
  load(file.path(data_dir, "state_inputs", paste0("pop_", fipscode, ".RData")), envir = e)
  load(file.path(data_dir, "state_inputs", paste0("le_", fipscode, ".RData")), envir = e)

  s <- list(
    M = list(initAC = e$m_M.initAC, cessAC = e$m_M.cessAC, NS.LE = e$m_M.NS.LE, census = e$df_M.census_data,
             mort = prepare_mortality(e$m_p_M.mortNS_AC, e$m_p_M.mortCS_AC, e$a_p_M.mortYSQ_AC),
             mort_ap = prepare_mort_ap(e$m_p_M.mortNS_AP, e$m_p_M.mortCS_AP, e$a_p_M.mortYSQ_AP, e$df_M.census_data)),
    F = list(initAC = e$m_F.initAC, cessAC = e$m_F.cessAC, NS.LE = e$m_F.NS.LE, census = e$df_F.census_data,
             mort = prepare_mortality(e$m_p_F.mortNS_AC, e$m_p_F.mortCS_AC, e$a_p_F.mortYSQ_AC),
             mort_ap = prepare_mort_ap(e$m_p_F.mortNS_AP, e$m_p_F.mortCS_AP, e$a_p_F.mortYSQ_AP, e$df_F.census_data)))
  for (g in c("M", "F")) {
    s[[g]]$base <- run_cohorts(s[[g]]$initAC, s[[g]]$cessAC, s[[g]], full = TRUE)
  }
  .cache$key <- key
  .cache$state <- s
  s
}


#===============================================================================
# cohort life-course smoking transitions, output prevalence
#===============================================================================

# Prevalence and population by age and cohort -> the list returned by the original generate_prevs().
# full = FALSE leaves out m_NSprevAP, which only the full runstates() outputs read.
assemble_prevs <- function(m_NSprevAC, m_CSprevAC, m_popAC, a_FSprevAC, full = TRUE) {
  cc <- as.character(startbc:endbc)
  colnames(m_NSprevAC) <- colnames(m_CSprevAC) <- colnames(m_popAC) <- cc

  to_ap <- function(m) { out <- matrix(0, 100, .nc + 100); out[.ap_lin] <- m; out }
  m_popAP <- to_ap(m_popAC)
  m_CSprevAP <- to_ap(m_CSprevAC)
  m_NSprevAP <- if (full) to_ap(m_NSprevAC)
  a_FSprevAP <- NULL
  if (!is.null(a_FSprevAC)) {
    a_FSprevAP <- array(0, c(100, .nc + 100, 40))
    a_FSprevAP[.ap_lin3] <- a_FSprevAC
    colnames(a_FSprevAP) <- as.character(startbc:endyear)
  }
  m_smokersAP <- m_CSprevAP * m_popAP
  colnames(m_popAP) <- colnames(m_CSprevAP) <- as.character(startbc:endyear)
  if (full) colnames(m_NSprevAP) <- as.character(startbc:endyear)

  list(m_output = NULL, m_NSprevAC = m_NSprevAC, m_CSprevAC = m_CSprevAC,
       a_FSprevAC = a_FSprevAC, m_popAC = m_popAC, m_NSprevAP = m_NSprevAP,
       m_CSprevAP = m_CSprevAP, a_FSprevAP = a_FSprevAP, m_popAP = m_popAP, m_smokersAP = m_smokersAP)
}

cohort_prevs <- function(m_init_AC, m_cess_AC, mort) {
  cc <- as.character(startbc:endbc)
  ip <- t(m_init_AC[, cc]); cp <- t(m_cess_AC[, cc])
  sNS <- mort$survNS_t; sCS <- mort$survCS_t; sY <- mort$survYSQ_t

  NS <- CS <- matrix(0, .nc, 100)
  Fm <- matrix(0, .nc * 40, 100)
  NS[, 1] <- v_stdbirths[1:.nc]
  for (i in 2:100) {
    ns <- NS[, i-1]; cs <- CS[, i-1]; p_i <- ip[, i-1]; p_c <- cp[, i-1]
    s_ns <- sNS[, i-1]; s_cs <- sCS[, i-1]; f <- Fm[, i-1]; s_y <- sY[, i-1]
    NS[, i] <- ns * (1 - p_i) * s_ns
    CS[, i] <- ns * p_i * s_ns + cs * (1 - p_c) * s_cs
    Fm[, i] <- c(cs * (p_c) * s_cs,
                 f[.fs_shift] * s_y[.fs_shift],
                 f[.fs_39] * s_y[.fs_39] + f[.fs_40] * s_y[.fs_40])
  }
  NS <- t(NS); CS <- t(CS)
  FS <- aperm(array(Fm, c(.nc, 40, 100)), c(3, 1, 2))
  m_popAC <- NS + CS + rowSums(FS, dims = 2)
  assemble_prevs(NS / m_popAC, CS / m_popAC, m_popAC, FS / as.vector(m_popAC))
}


#===============================================================================
# Smoking-attributable deaths (SADs) and years of life lost (YLL)
#===============================================================================

# The original AP-to-AC loop assigns m_AC[i, j - i + 1] for every j, so for j < i - 1 the column
# index is negative and fills every other column of row i. The cells it never assigns afterwards
# (row i >= 3, columns past 193 - i + 1) therefore hold m_AP[i, i - 2]; reproduced here exactly.
ap_to_ac <- function(m_AP) {
  m_AC <- matrix(NA_real_, 100, .nc)
  m_AC[.ac_ok_dst] <- m_AP[.ac_ok_src]
  m_AC[.ac_tail_dst] <- m_AP[.ac_tail_src]
  m_AC
}

sad_ap <- function(l_prev_outputs, mort_ap) {
  m_popdist <- mort_ap$popdist
  m_SAD_AP <- m_popdist * (l_prev_outputs$m_CSprevAP[, v_calyears] * mort_ap$dCS)
  for (j in 1:40) {
    m_SAD_AP <- m_SAD_AP + m_popdist * l_prev_outputs$a_FSprevAP[, v_calyears, j] * mort_ap$dYSQ[, , j]
  }
  m_SAD_AP
}

# full = FALSE returns df_SAD_AP and df_YLL_AP as NULL; the vectors and AC matrices are unchanged.
mort_outcomes <- function(m_SAD_AP, m_NS.LE, full = TRUE) {
  if (full) {
    df_SAD_AP <- as.data.frame(m_SAD_AP)
    v_SADyear <- colSums(df_SAD_AP)
    df_YLL_AP <- df_SAD_AP * m_NS.LE
    v_YLLyear <- colSums(df_YLL_AP)
    m_SAD_AP <- as.matrix(df_SAD_AP)
    m_YLL_AP <- as.matrix(df_YLL_AP)
  } else {
    df_SAD_AP <- df_YLL_AP <- NULL
    v_SADyear <- colSums(m_SAD_AP)
    m_YLL_AP <- m_SAD_AP * m_NS.LE
    v_YLLyear <- colSums(m_YLL_AP)
  }

  m_SAD_AC <- ap_to_ac(m_SAD_AP)
  m_YLL_AC <- ap_to_ac(m_YLL_AP)
  colnames(m_SAD_AC) <- colnames(m_YLL_AC) <- as.character(startbc:endbc)
  rownames(m_SAD_AC) <- rownames(m_YLL_AC) <- 0:99

  list(df_SAD_AP = df_SAD_AP, v_SADyear = v_SADyear,
       df_YLL_AP = df_YLL_AP, v_YLLyear = v_YLLyear,
       m_SAD_AC = m_SAD_AC, m_YLL_AC = m_YLL_AC)
}

# Cohort model plus mortality outcomes for one sex. g holds that sex's cached state inputs.
# full = FALSE leaves out the former-smoker prevalence arrays, which only feed the SAD sum, and the
# data-frame copies of SAD and YLL.
run_cohorts <- function(m_init_AC, m_cess_AC, g, full = TRUE) {
  if (use_kernel()) {
    cc <- as.character(startbc:endbc)
    k <- .kernel$cohort_kernel(m_init_AC[, cc], m_cess_AC[, cc], g$mort$survNS_ac, g$mort$survCS_ac,
                               g$mort$survYSQ_jt, v_stdbirths[1:.nc], g$mort_ap$popdist,
                               g$mort_ap$dCS, g$mort_ap$dYSQ_jt, full)
    prev <- assemble_prevs(k$NSprev, k$CSprev, k$pop, k$FSprev, full)
    m_SAD_AP <- k$SAD
    dimnames(m_SAD_AP) <- dimnames(g$mort_ap$popdist)
  } else {
    prev <- cohort_prevs(m_init_AC, m_cess_AC, g$mort)
    m_SAD_AP <- sad_ap(prev, g$mort_ap)
    if (!full) prev["a_FSprevAC"] <- prev["a_FSprevAP"] <- list(NULL)
  }
  list(prev = prev, mort = mort_outcomes(m_SAD_AP, g$NS.LE, full))
}


#===============================================================================
# Run one state under one set of policy effects
#===============================================================================

# full = FALSE returns df_mort.outputs and df_CSprevs.by.state as NULL and leaves the
# former-smoker prevalence arrays out of the policy prevalence objects; the TCP tool outputs
# (generate_TCPoutput) are unchanged.
runstates <- function(fipscode, m.initiation.effect, m.cessation.effect, data_dir = "data",
                      full = TRUE) {
  s <- get_state(fipscode, data_dir)

  policy_AC <- function(m_AC, effect) {
    m <- matrix(m_AC[.ac_lin] * effect[.ap_lin], 100, .nc)
    colnames(m) <- as.character(startbc:endbc)
    m
  }
  male <- run_cohorts(policy_AC(s$M$initAC, m.initiation.effect),
                      policy_AC(s$M$cessAC, m.cessation.effect), s$M, full)
  female <- run_cohorts(policy_AC(s$F$initAC, m.initiation.effect),
                        policy_AC(s$F$cessAC, m.cessation.effect), s$F, full)
  l_M.policy.prev <- male$prev; l_M.policy.mort <- male$mort
  l_F.policy.prev <- female$prev; l_F.policy.mort <- female$mort
  l_M.base.prev <- s$M$base$prev; l_F.base.prev <- s$F$base$prev
  l_M.base.mort <- s$M$base$mort; l_F.base.mort <- s$F$base$mort
  df_M.census_data <- s$M$census; df_F.census_data <- s$F$census

  l_prev_out <- list(
    state = fipscode,
    m_M_smokers = l_M.policy.prev$m_smokersAP,
    m_F_smokers = l_F.policy.prev$m_smokersAP,
    m_M_popAP = l_M.policy.prev$m_popAP,
    m_F_popAP = l_F.policy.prev$m_popAP,
    male_prev_obj        = l_M.policy.prev,
    female_prev_obj      = l_F.policy.prev,
    male_base_prev_obj   = l_M.base.prev,
    female_base_prev_obj = l_F.base.prev,
    male_mort_obj        = l_M.policy.mort,
    female_mort_obj      = l_F.policy.mort,
    male_base_mort_obj   = l_M.base.mort,
    female_base_mort_obj = l_F.base.mort,
    startbc   = startbc,
    endbc     = endbc,
    endyear   = endyear
  )
  if (!full) {
    return(list(df_mort.outputs = NULL, l_prev_out = l_prev_out, df_CSprevs.by.state = NULL))
  }

  #------------------- format prev for outputting -----------------------------------
  df_CSprevs.by.state <- NULL
  v_minage <- c(18, 18, 25, 45, 65)
  v_maxage <- c(99, 24, 44, 64, 99)
  m_M.CSprev <- l_M.policy.prev$m_CSprevAP
  m_F.CSprev <- l_F.policy.prev$m_CSprevAP
  m.M.pop_AP <- as.matrix(cbind(df_M.census_data, rep(df_M.census_data[cohyears], 100)))
  m.F.pop_AP <- as.matrix(cbind(df_F.census_data, rep(df_F.census_data[cohyears], 100)))
  colnames(m.M.pop_AP) <- colnames(m.F.pop_AP) <- as.character(startbc:endyear)

  for (i in 1:5) {
    rows <- (v_minage[i] + 1):(v_maxage[i] + 1)
    v_M.num <- colSums(m.M.pop_AP[rows, ] * m_M.CSprev[rows, ])
    v_F.num <- colSums(m.F.pop_AP[rows, ] * m_F.CSprev[rows, ])
    v_M.den <- colSums(m.M.pop_AP[rows, ])
    v_F.den <- colSums(m.F.pop_AP[rows, ])
    v_M.prev.minmax <- v_M.num / v_M.den
    v_F.prev.minmax <- v_F.num / v_F.den
    v_B.prev.minmax <- (v_M.num + v_F.num) / (v_M.den + v_F.den)

    df_CSprevbystate_temp <- as.data.frame(rbind(
      cbind(v_M.prev.minmax[1:cohyears], "Men"),
      cbind(v_F.prev.minmax[1:cohyears], "Women"),
      cbind(v_B.prev.minmax[1:cohyears], "Both")
    ))
    colnames(df_CSprevbystate_temp) <- c("prev", "gender")
    df_CSprevbystate_temp$age <- paste0(v_minage[i], ".", v_maxage[i])
    df_CSprevbystate_temp$year <- rep(names(v_M.prev.minmax[1:cohyears]), 3)
    df_CSprevbystate_temp$state <- fipscode
    df_CSprevbystate_temp$abbr <- fips_to_abbr(fipscode)
    df_CSprevs.by.state <- rbind(df_CSprevs.by.state, df_CSprevbystate_temp)
  }
  df_CSprevs.by.state$prev <- as.numeric(df_CSprevs.by.state$prev)
  df_CSprevs.by.state$year <- as.numeric(df_CSprevs.by.state$year)

  #---------------format policy YLL and SADs for outputting --------------------
  v_M.YLLyear <- l_M.policy.mort$v_YLLyear
  v_F.YLLyear <- l_F.policy.mort$v_YLLyear
  v_B.YLLyear <- v_M.YLLyear + v_F.YLLyear
  v_M.SADyear <- l_M.policy.mort$v_SADyear
  v_F.SADyear <- l_F.policy.mort$v_SADyear
  v_B.SADyear <- v_M.SADyear + v_F.SADyear

  v_M.LYGyear <- colSums(l_M.base.mort$df_YLL_AP - l_M.policy.mort$df_YLL_AP)
  v_F.LYGyear <- colSums(l_F.base.mort$df_YLL_AP - l_F.policy.mort$df_YLL_AP)
  v_B.LYGyear <- v_M.LYGyear + v_F.LYGyear

  v_M.SADs_averted_year <- l_M.base.mort$v_SADyear - l_M.policy.mort$v_SADyear
  v_F.SADs_averted_year <- l_F.base.mort$v_SADyear - l_F.policy.mort$v_SADyear
  v_B.SADs_averted_year <- v_M.SADs_averted_year + v_F.SADs_averted_year

  mortout <- function(YLL, SAD, LYG, SADav, label) {
    cbind(YLL, cumsum(YLL), SAD, cumsum(SAD), LYG, cumsum(LYG), SADav, cumsum(SADav), startbc:endbc, label)
  }
  df_mort.outputs <- as.data.frame(rbind(
    mortout(v_M.YLLyear, v_M.SADyear, v_M.LYGyear, v_M.SADs_averted_year, "Male"),
    mortout(v_F.YLLyear, v_F.SADyear, v_F.LYGyear, v_F.SADs_averted_year, "Female"),
    mortout(v_B.YLLyear, v_B.SADyear, v_B.LYGyear, v_B.SADs_averted_year, "Both")))
  colnames(df_mort.outputs) <- c("YLL", "YLLcum", "SADs", "SADcum", "LYG", "LYGcum",
                                 "SADsAverted", "SADsAvertedcum", "year", "gender")
  df_mort.outputs$state <- fipscode
  df_mort.outputs$abbr <- fips_to_abbr(fipscode)
  num_cols <- c("YLL", "YLLcum", "SADs", "SADcum", "LYG", "LYGcum", "SADsAverted", "SADsAvertedcum", "year")
  df_mort.outputs[num_cols] <- lapply(df_mort.outputs[num_cols], as.numeric)

  list(df_mort.outputs = df_mort.outputs,
       l_prev_out = l_prev_out,
       df_CSprevs.by.state = df_CSprevs.by.state)
}


#===============================================================================
# TCP Tool output tables
#===============================================================================

tcp_results <- function(l_prev_out, policyyear, coh_keep) {
  mb <- l_prev_out$male_base_prev_obj; fb <- l_prev_out$female_base_prev_obj
  mp <- l_prev_out$male_prev_obj;      fp <- l_prev_out$female_prev_obj
  yrs <- policyyear:l_prev_out$endyear
  col_of <- match(as.character(yrs), colnames(mb$m_popAP))
  if (anyNA(col_of)) stop("Year ", yrs[is.na(col_of)][1], " not found in AP matrices.")

  cohs <- as.integer(coh_keep)
  year <- rep(yrs, each = length(cohs))
  cohort <- rep(cohs, times = length(yrs))
  age <- year - cohort
  keep <- age >= 0 & age <= 99
  year <- year[keep]; cohort <- cohort[keep]; age <- age[keep]
  ix <- cbind(age + 1, col_of[match(year, yrs)])
  df_cohort <- data.frame(
    year = year, age = as.character(age), cohort = as.character(cohort),
    males_baseline = mb$m_smokersAP[ix] / mb$m_popAP[ix],
    females_baseline = fb$m_smokersAP[ix] / fb$m_popAP[ix],
    males_policy = mp$m_smokersAP[ix] / mp$m_popAP[ix],
    females_policy = fp$m_smokersAP[ix] / fp$m_popAP[ix],
    both_baseline = (mb$m_smokersAP[ix] + fb$m_smokersAP[ix]) / (mb$m_popAP[ix] + fb$m_popAP[ix]),
    both_policy = (mp$m_smokersAP[ix] + fp$m_smokersAP[ix]) / (mp$m_popAP[ix] + fp$m_popAP[ix]),
    policy_year = policyyear,
    stringsAsFactors = FALSE)

  agegroups <- c("12-17", "18-24", "25-44", "45-64", "65p", "18-99")
  agel <- c(12, 18, 25, 45, 65, 18)
  ageu <- c(17, 24, 44, 64, 99, 99)
  safe_rate <- function(n, d) ifelse(d > 0, n / d, 0)
  cs <- function(m, rows) colSums(m[rows, col_of, drop = FALSE])
  df_groups <- do.call(rbind, lapply(seq_along(agegroups), function(kag) {
    rows <- (agel[kag]:ageu[kag]) + 1
    Ms_b <- cs(mb$m_smokersAP, rows); Mp_b <- cs(mb$m_popAP, rows)
    Fs_b <- cs(fb$m_smokersAP, rows); Fp_b <- cs(fb$m_popAP, rows)
    Ms_p <- cs(mp$m_smokersAP, rows); Mp_p <- cs(mp$m_popAP, rows)
    Fs_p <- cs(fp$m_smokersAP, rows); Fp_p <- cs(fp$m_popAP, rows)
    data.frame(year = yrs, age = agegroups[kag], cohort = "ALL",
               males_baseline = safe_rate(Ms_b, Mp_b), females_baseline = safe_rate(Fs_b, Fp_b),
               males_policy = safe_rate(Ms_p, Mp_p), females_policy = safe_rate(Fs_p, Fp_p),
               both_baseline = safe_rate(Ms_b + Fs_b, Mp_b + Fp_b),
               both_policy = safe_rate(Ms_p + Fs_p, Mp_p + Fp_p),
               policy_year = policyyear, stringsAsFactors = FALSE)
  }))
  df_groups <- df_groups[order(match(df_groups$year, yrs), match(df_groups$age, agegroups)), ]

  out <- rbind(df_cohort, df_groups)
  rownames(out) <- NULL
  out[order(out$cohort, out$age, out$year), ]
}

generate_TCPoutput <- function(l_prev_out, fipscode, policyyear,
                               s_cohorts = seq(1970, 2030, by = 10),
                               out_root = "tcp_tool",
                               policy_name = NULL,
                               scenario_tag = "scenario",
                               write_files = TRUE) {
  if (is.null(policy_name) || !nzchar(policy_name)) {
    return(invisible(NULL))
  }
  abbr <- fips_to_abbr(fipscode)

  coh_keep <- as.character(s_cohorts)
  coh_keep <- coh_keep[coh_keep %in% colnames(l_prev_out$male_base_mort_obj$m_SAD_AC)]
  if (length(coh_keep) == 0) stop("None of s_cohorts exist in the AC matrices.")

  cum_by_age <- function(m_base_AC, m_pol_AC) {
    m_diff <- m_base_AC[, coh_keep, drop = FALSE] - m_pol_AC[, coh_keep, drop = FALSE]
    m_diff[is.na(m_diff)] <- 0
    apply(m_diff, 2, cumsum)
  }
  make_cohort_long <- function(M_mat, F_mat, value_M, value_F) {
    age <- rep(0:99, times = length(coh_keep))
    cohort <- rep(coh_keep, each = 100)
    out <- data.frame(year = as.integer(cohort) + age, cohort = cohort, stringsAsFactors = FALSE)
    out[[value_M]] <- as.numeric(M_mat)
    out[[value_F]] <- as.numeric(F_mat)
    out[out$year >= policyyear, , drop = FALSE]
  }
  make_all <- function(v_base_M, v_pol_M, v_base_F, v_pol_F, col_M, col_F) {
    yrs <- as.integer(names(v_base_M))
    if (length(yrs) == 0) stop("Annual vectors must have names equal to calendar years.")
    ord <- order(yrs)
    df <- data.frame(year = yrs[ord], cohort = "ALL", stringsAsFactors = FALSE)
    df[[col_M]] <- cumsum(as.numeric(v_base_M)[ord] - as.numeric(v_pol_M)[ord])
    df[[col_F]] <- cumsum(as.numeric(v_base_F)[ord] - as.numeric(v_pol_F)[ord])
    df[df$year >= policyyear, , drop = FALSE]
  }

  p <- l_prev_out
  deathsfile <- rbind(
    make_cohort_long(cum_by_age(p$male_base_mort_obj$m_SAD_AC, p$male_mort_obj$m_SAD_AC),
                     cum_by_age(p$female_base_mort_obj$m_SAD_AC, p$female_mort_obj$m_SAD_AC),
                     "deaths_avoided_males", "deaths_avoided_females"),
    make_all(p$male_base_mort_obj$v_SADyear, p$male_mort_obj$v_SADyear,
             p$female_base_mort_obj$v_SADyear, p$female_mort_obj$v_SADyear,
             "deaths_avoided_males", "deaths_avoided_females"))
  lygfile <- rbind(
    make_cohort_long(cum_by_age(p$male_base_mort_obj$m_YLL_AC, p$male_mort_obj$m_YLL_AC),
                     cum_by_age(p$female_base_mort_obj$m_YLL_AC, p$female_mort_obj$m_YLL_AC),
                     "cLYG_males", "cLYG_females"),
    make_all(p$male_base_mort_obj$v_YLLyear, p$male_mort_obj$v_YLLyear,
             p$female_base_mort_obj$v_YLLyear, p$female_mort_obj$v_YLLyear,
             "cLYG_males", "cLYG_females"))
  deathsfile$policy_year <- policyyear
  lygfile$policy_year <- policyyear
  resultsfile <- tcp_results(l_prev_out, policyyear, coh_keep)

  if (write_files) {
    policy_part <- if (!is.null(policy_name) && nzchar(policy_name)) policy_name else NULL
    dirs <- lapply(c("deaths", "lyg", "results"), function(d) do.call(file.path, as.list(c(out_root, abbr, policy_part, d))))
    for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(deathsfile,  file.path(dirs[[1]], paste0("deaths_",  scenario_tag, ".csv")))
    readr::write_csv(lygfile,     file.path(dirs[[2]], paste0("lyg_",     scenario_tag, ".csv")))
    readr::write_csv(resultsfile, file.path(dirs[[3]], paste0("results_", scenario_tag, ".csv")))
  }
  invisible(list(deaths = deathsfile, lyg = lygfile, results = resultsfile))
}


# The tables generate_TCPoutput(write_files = FALSE) returns, for many scenarios of one state.
# Everything that is the same for every scenario of the state (row layout and order, the baseline
# columns and baseline sums) is computed on the first call for each policy year and kept in
# `memo`, an environment the caller keeps per state; later calls compute only the policy values.
# Returns, for deaths, lyg and results, list(labels = <data frame of the non-value columns>,
# values = <numeric matrix of the value columns>); tcp_table_df() joins them into the data frame
# generate_TCPoutput() returns.
tcp_tables <- function(l_prev_out, policyyear, memo, s_cohorts = seq(1970, 2030, by = 10)) {
  p <- l_prev_out
  key <- as.character(policyyear)
  m <- memo[[key]]
  if (is.null(m)) {
    m <- tcp_layout(p, policyyear, s_cohorts)
    assign(key, m, envir = memo)
  }

  # deaths and lyg: cumulative baseline - policy by age within each kept cohort, then by year
  cum_long <- function(b_AC, p_AC) {
    d <- b_AC - p_AC[, m$coh_keep, drop = FALSE]
    d[is.na(d)] <- 0
    as.numeric(apply(d, 2, cumsum))[m$long_keep]
  }
  cum_all <- function(b_v, p_v) cumsum(b_v - as.numeric(p_v)[m$all_ord])[m$all_keep]
  mort_values <- function(what, cols) {
    v <- cbind(c(cum_long(m$base[[paste0("M_", what, "_AC")]], p$male_mort_obj[[paste0("m_", what, "_AC")]]),
                 cum_all(m$base[[paste0("M_", what, "year")]], p$male_mort_obj[[paste0("v_", what, "year")]])),
               c(cum_long(m$base[[paste0("F_", what, "_AC")]], p$female_mort_obj[[paste0("m_", what, "_AC")]]),
                 cum_all(m$base[[paste0("F_", what, "year")]], p$female_mort_obj[[paste0("v_", what, "year")]])))
    colnames(v) <- cols
    v
  }

  # results: policy prevalence for the cohort rows and the age-group rows
  mp <- p$male_prev_obj; fp <- p$female_prev_obj
  Ms <- mp$m_smokersAP[m$cell]; Mp <- mp$m_popAP[m$cell]
  Fs <- fp$m_smokersAP[m$cell]; Fp <- fp$m_popAP[m$cell]
  safe_rate <- function(n, d) ifelse(d > 0, n / d, 0)
  cs <- function(x, rows) colSums(x[rows, m$col_of, drop = FALSE])
  grp <- lapply(m$group_rows, function(rows) {
    Ms_p <- cs(mp$m_smokersAP, rows); Mp_p <- cs(mp$m_popAP, rows)
    Fs_p <- cs(fp$m_smokersAP, rows); Fp_p <- cs(fp$m_popAP, rows)
    cbind(safe_rate(Ms_p, Mp_p), safe_rate(Fs_p, Fp_p), safe_rate(Ms_p + Fs_p, Mp_p + Fp_p))
  })
  pol <- rbind(cbind(Ms / Mp, Fs / Fp, (Ms + Fs) / (Mp + Fp)),
               do.call(rbind, grp)[m$group_ord, , drop = FALSE])[m$results_ord, , drop = FALSE]
  res <- cbind(m$base$results[, 1], m$base$results[, 2], pol[, 1], pol[, 2], m$base$results[, 3], pol[, 3])
  colnames(res) <- c("males_baseline", "females_baseline", "males_policy", "females_policy",
                     "both_baseline", "both_policy")

  list(deaths = list(labels = m$labels$deaths, values = mort_values("SAD", c("deaths_avoided_males", "deaths_avoided_females"))),
       lyg = list(labels = m$labels$lyg, values = mort_values("YLL", c("cLYG_males", "cLYG_females"))),
       results = list(labels = m$labels$results, values = res))
}

# The scenario-independent part of tcp_tables() for one state and policy year.
tcp_layout <- function(p, policyyear, s_cohorts) {
  coh_keep <- as.character(s_cohorts)
  coh_keep <- coh_keep[coh_keep %in% colnames(p$male_base_mort_obj$m_SAD_AC)]
  if (length(coh_keep) == 0) stop("None of s_cohorts exist in the AC matrices.")

  # deaths / lyg rows: each kept cohort by age from the policy year, then ALL by year
  long_year <- as.integer(rep(coh_keep, each = 100)) + rep(0:99, times = length(coh_keep))
  long_keep <- long_year >= policyyear
  all_years <- as.integer(names(p$male_base_mort_obj$v_SADyear))
  all_ord <- order(all_years)
  all_keep <- all_years[all_ord] >= policyyear
  mort_labels <- data.frame(year = c(long_year[long_keep], all_years[all_ord][all_keep]),
                            cohort = c(rep(coh_keep, each = 100)[long_keep], rep("ALL", sum(all_keep))),
                            policy_year = policyyear, stringsAsFactors = FALSE)
  base <- list()
  for (g in c("M", "F")) {
    b <- p[[if (g == "M") "male_base_mort_obj" else "female_base_mort_obj"]]
    for (what in c("SAD", "YLL")) {
      base[[paste0(g, "_", what, "_AC")]] <- b[[paste0("m_", what, "_AC")]][, coh_keep, drop = FALSE]
      base[[paste0(g, "_", what, "year")]] <- as.numeric(b[[paste0("v_", what, "year")]])[all_ord]
    }
  }

  # results rows: (year, cohort) cells, then age groups by year, sorted as tcp_results() sorts
  mb <- p$male_base_prev_obj; fb <- p$female_base_prev_obj
  yrs <- policyyear:p$endyear
  col_of <- match(as.character(yrs), colnames(mb$m_popAP))
  if (anyNA(col_of)) stop("Year ", yrs[is.na(col_of)][1], " not found in AP matrices.")
  cohs <- as.integer(coh_keep)
  year <- rep(yrs, each = length(cohs))
  cohort <- rep(cohs, times = length(yrs))
  age <- year - cohort
  keep <- age >= 0 & age <= 99
  year <- year[keep]; cohort <- cohort[keep]; age <- age[keep]
  cell <- (col_of[match(year, yrs)] - 1) * nrow(mb$m_popAP) + age + 1
  agegroups <- c("12-17", "18-24", "25-44", "45-64", "65p", "18-99")
  agel <- c(12, 18, 25, 45, 65, 18)
  ageu <- c(17, 24, 44, 64, 99, 99)
  group_rows <- lapply(seq_along(agegroups), function(k) (agel[k]:ageu[k]) + 1)
  g_year <- rep(yrs, times = length(agegroups))
  g_age <- rep(agegroups, each = length(yrs))
  group_ord <- order(match(g_year, yrs), match(g_age, agegroups))
  res_labels <- data.frame(year = c(year, g_year[group_ord]), age = c(as.character(age), g_age[group_ord]),
                           cohort = c(as.character(cohort), rep("ALL", length(g_year))),
                           policy_year = policyyear, stringsAsFactors = FALSE)
  results_ord <- order(res_labels$cohort, res_labels$age, res_labels$year)
  res_labels <- res_labels[results_ord, ]
  rownames(res_labels) <- NULL

  Ms <- mb$m_smokersAP[cell]; Mp <- mb$m_popAP[cell]
  Fs <- fb$m_smokersAP[cell]; Fp <- fb$m_popAP[cell]
  safe_rate <- function(n, d) ifelse(d > 0, n / d, 0)
  cs <- function(x, rows) colSums(x[rows, col_of, drop = FALSE])
  grp <- lapply(group_rows, function(rows) {
    Ms_b <- cs(mb$m_smokersAP, rows); Mp_b <- cs(mb$m_popAP, rows)
    Fs_b <- cs(fb$m_smokersAP, rows); Fp_b <- cs(fb$m_popAP, rows)
    cbind(safe_rate(Ms_b, Mp_b), safe_rate(Fs_b, Fp_b), safe_rate(Ms_b + Fs_b, Mp_b + Fp_b))
  })
  base$results <- rbind(cbind(Ms / Mp, Fs / Fp, (Ms + Fs) / (Mp + Fp)),
                        do.call(rbind, grp)[group_ord, , drop = FALSE])[results_ord, , drop = FALSE]

  list(coh_keep = coh_keep, long_keep = long_keep, all_ord = all_ord, all_keep = all_keep,
       col_of = col_of, cell = cell, group_rows = group_rows, group_ord = group_ord,
       results_ord = results_ord, base = base,
       labels = list(deaths = mort_labels, lyg = mort_labels, results = res_labels))
}

# A tcp_tables() table as the data frame generate_TCPoutput() returns (policy_year last).
tcp_table_df <- function(tab) {
  lab <- tab$labels
  df <- lab[setdiff(names(lab), "policy_year")]
  for (cl in colnames(tab$values)) df[[cl]] <- tab$values[, cl]
  df$policy_year <- lab$policy_year
  df
}


#===============================================================================
# Cigarette Tax
#===============================================================================

tax_effectCalculation <- function(initprice, tax,
                                  startbc = 1908, endyear = 2200, policyYear,
                                  inidecay = 0.0, cesdecay = 0.2, iniagemod = 1, cesagemod = 1) {
  ages <- 0:99
  periods <- startbc:endyear
  cess_elasticities <- c(rep(0, 10), rep(-2.00, 90))
  init_elasticities <- c(rep(0, 10), rep(-0.4, 8), rep(-0.3, 7), rep(-0.2, 20), rep(0, 55))

  newprice <- initprice + tax
  pricechange <- (newprice - initprice) / ((newprice + initprice) / 2)

  m.initiation.effect <- matrix(1, nrow = length(ages), ncol = length(periods))
  m.cessation.effect  <- matrix(1, nrow = length(ages), ncol = length(periods))
  colnames(m.initiation.effect) <- colnames(m.cessation.effect) <- as.character(periods)

  post <- which(periods >= policyYear)
  timeSincePolicy <- periods[post] - policyYear
  m.initiation.effect[, post] <- 1 - (-pricechange * init_elasticities)
  m.cessation.effect[, post] <- 1 + outer(-pricechange * cess_elasticities, (1 - cesdecay)^timeSincePolicy)

  list(m.initiation.effect = m.initiation.effect,
       m.cessation.effect  = m.cessation.effect)
}
