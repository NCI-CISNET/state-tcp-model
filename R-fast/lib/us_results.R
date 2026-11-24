# US results for the TCP tool, computed the way update_us_results.py does it: each state's results
# values (rounded to the 6 significant digits the state files hold) weighted by that state's
# standardized model population, summed over states and divided by the summed weights. Plain R
# arithmetic, so values can differ from the script's in the last bits of a double, which does not
# change the 6-digit files except, rarely, by one in the last digit.

us_age_bands <- list("12-17" = c(12, 17), "18-24" = c(18, 24), "25-44" = c(25, 44),
                     "45-64" = c(45, 64), "18-99" = c(18, 99), "65p" = c(65, Inf))
us_age_group_sort <- c("12-17" = 0, "18-24" = 1, "25-44" = 2, "45-64" = 3, "65p" = 4, "18-99" = 5)

# Population weights by state and sex: `exact` keyed by "cohort_year year age" and `band` keyed
# by "year age_label" (sums over the band's ages). Cached next to the compiled code, keyed on the
# file's path, size and modification time.
load_us_population <- function(path, cache_dir = tools::R_user_dir("tcp_fast", "cache")) {
  path <- normalizePath(path)
  info <- file.info(path)
  key <- paste("v4", path, info$size, as.numeric(info$mtime))
  cache_file <- file.path(cache_dir, "us_population.rds")
  if (file.exists(cache_file)) {
    cached <- readRDS(cache_file)
    if (identical(cached$key, key)) return(cached$pop)
  }
  df <- vroom::vroom(path, delim = ",", col_types = "cccciiic", progress = FALSE, num_threads = 1)
  w <- as.numeric(df$model_population)
  sex <- tolower(df$sex)
  exact_key <- paste(df$cohort_year, df$year, df$age)
  rows <- split(seq_len(nrow(df)), list(df$abbr, sex), sep = "|")
  pop <- list()
  for (st in unique(df$abbr)) {
    pop[[st]] <- list()
    for (sx in c("male", "female")) {
      i <- rows[[paste(st, sx, sep = "|")]]
      pop[[st]][[sx]] <- list(exact_key = exact_key[i], exact = w[i])
    }
  }
  # Band sums by (state, sex, year); missing populations count as 0, as pandas sums them.
  bands <- do.call(rbind, lapply(names(us_age_bands), function(label) {
    r <- us_age_bands[[label]]
    i <- which(df$age >= r[1] & df$age <= r[2])
    grp <- paste(df$abbr[i], sex[i], df$year[i])
    u <- unique(grp)
    first <- i[match(u, grp)]
    data.frame(abbr = df$abbr[first], sex = sex[first], key = paste(df$year[first], label),
               sum = as.numeric(rowsum(w[i], match(grp, u), na.rm = TRUE)))
  }))
  brows <- split(seq_len(nrow(bands)), list(bands$abbr, bands$sex), sep = "|")
  for (st in names(pop)) {
    for (sx in c("male", "female")) {
      i <- brows[[paste(st, sx, sep = "|")]]
      pop[[st]][[sx]][["band_key"]] <- bands$key[i]
      pop[[st]][[sx]][["band"]] <- bands$sum[i]
    }
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(key = key, pop = pop), cache_file)
  pop
}

# Male and female population weight for each row of a state's results table (NA where the script's
# lookup finds none).
us_row_weights <- function(pop_st, res) {
  age <- as.character(res$age); cohort <- as.character(res$cohort)
  exact <- grepl("^[0-9]+$", age) & grepl("^[0-9]+$", cohort)
  ekey <- paste(cohort, res$year, age)
  bkey <- paste(res$year, age)
  lapply(c(male = "male", female = "female"), function(sx) {
    p <- pop_st[[sx]]
    out <- rep(NA_real_, nrow(res))
    if (is.null(p)) return(out)
    out[exact] <- p[["exact"]][match(ekey[exact], p[["exact_key"]])]
    if (length(p[["band"]])) out[!exact] <- p[["band"]][match(bkey[!exact], p[["band_key"]])]
    out
  })
}

# One state's contribution for each results row: weighted values and weights for male/female
# baseline/policy. `vals` holds the values the script reads from the state's results file.
us_terms <- function(vals, wts) {
  term <- function(v, p) {
    ok <- !is.na(v) & !is.na(p) & p > 0
    num <- v * p; num[!ok] <- 0
    wt <- p; wt[!ok] <- 0
    cbind(num, wt)
  }
  m <- cbind(term(vals$males_baseline, wts$male), term(vals$females_baseline, wts$female),
             term(vals$males_policy, wts$male), term(vals$females_policy, wts$female))
  colnames(m) <- c("mb_num", "mb_w", "fb_num", "fb_w", "mp_num", "mp_w", "fp_num", "fp_w")
  m
}

# Sum of the states' term matrices.
us_fold <- function(term_list) Reduce(`+`, term_list)

# The US results table from the summed terms, rows in the script's sort order. `labels` is the
# year/age/cohort/policy_year layout shared by every state's results table.
us_results_table <- function(s, labels) {
  ratio <- function(num, den) ifelse(den > 0, num, 0) / ifelse(den > 0, den, 1)
  bw <- s[, "mb_w"] + s[, "fb_w"]
  pw <- s[, "mp_w"] + s[, "fp_w"]
  out <- data.frame(year = labels$year, age = as.character(labels$age), cohort = as.character(labels$cohort),
                    males_baseline = ratio(s[, "mb_num"], s[, "mb_w"]),
                    females_baseline = ratio(s[, "fb_num"], s[, "fb_w"]),
                    males_policy = ratio(s[, "mp_num"], s[, "mp_w"]),
                    females_policy = ratio(s[, "fp_num"], s[, "fp_w"]),
                    both_baseline = ratio(s[, "mb_num"] + s[, "fb_num"], bw),
                    both_policy = ratio(s[, "mp_num"] + s[, "fp_num"], pw),
                    policy_year = labels$policy_year, stringsAsFactors = FALSE)
  for (cl in c("males_baseline", "females_baseline", "males_policy", "females_policy", "both_baseline", "both_policy")) {
    out[[cl]][is.na(out[[cl]])] <- 0
  }
  cohort_num <- suppressWarnings(as.numeric(ifelse(grepl("^[0-9]+$", out$cohort), out$cohort, NA)))
  age_num <- suppressWarnings(as.numeric(ifelse(grepl("^[0-9]+$", out$age), out$age, NA)))
  age_group <- unname(us_age_group_sort[out$age]); age_group[is.na(age_group)] <- 999
  o <- order(out$policy_year, is.na(cohort_num), ifelse(is.na(cohort_num), 999999, cohort_num),
             is.na(age_num), age_group, ifelse(is.na(age_num), 999999, age_num), out$year, method = "radix")
  out <- out[o, ]
  rownames(out) <- NULL
  out
}
