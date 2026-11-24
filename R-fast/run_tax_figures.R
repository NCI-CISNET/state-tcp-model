# Tax policy figures for every state, with the fast model in lib/tcp_fast.R.
#
# Runs the analysis behind R/04_visualization_policy_specific.R's tax figures, as
# 03_main_analysis.R on main defined it: for each state, the no-policy baseline and tax increases
# of +$1 ... +$5 (figures.tax_increases) starting in figures.policy_year, each with full
# runstates() outputs. Saves, in figures.out_dir:
#   baseline_results.rds  list keyed "<fips>" with df_mort.outputs and df_CSprevs.by.state
#   tax_out_list.rds      list keyed "<fips>_<d>dollar" with df_mort.outputs, df_CSprevs.by.state,
#                         initprice and tax_raise
#   Both lists also have US entries ("US", "US_<d>dollar"): all states together, each with its own
#   tax increase. Deaths and life-years are summed over states; prevalence is weighted by each
#   state's census population, as runstates() weights ages within a state.
#   figs/tax/<fips>_<State>_{LYG,DeathAvoided,Prev,PrevRed}_tax.<svg|jpg>   (figures.formats), and
#                         US_United_States_*; each state's prevalence figures show the US in grey
#   index.html            every figure on one page, by state (needs svg in figures.formats)
# figures.plots = off writes only the .rds files.
# The figures come from R/tax_plots.R, the same functions R/04 uses. The lists keep only the parts
# of each runstates() result the figure scripts read; the original also kept l_prev_out (the
# full prevalence arrays, ~50 MB per scenario), which nothing downstream uses.
#
# usage: Rscript R-fast/run_tax_figures.R [--config=<config.yml>] [--<setting>=<value> ...]
#   Settings: data-dir, prices (inputs); out-dir, policy-year, tax-increases, states, plots,
#   formats, html, workers (figures). See config.example.yml.
fast_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])))
source(file.path(fast_dir, "lib", "packages.R"))
# dplyr, ggplot2, scales and svglite are for the figures (R/tax_plots.R).
use_packages(c("Rcpp", "usmap", "yaml", "dplyr", "ggplot2", "scales", "svglite"))
suppressMessages({ library(parallel) })
repo <- dirname(fast_dir)
source(file.path(fast_dir, "lib", "config.R"))
source(file.path(fast_dir, "lib", "tcp_fast.R"))

settings <- read_settings(file.path(fast_dir, "config.yml"))
data_dir <- normalizePath(settings$get("inputs", "data_dir", type = "path", required = TRUE))
prices_csv <- settings$get("inputs", "prices", type = "path", required = TRUE)
out_dir <- settings$get("figures", "out_dir", type = "path", required = TRUE)
policy_year <- settings$get("figures", "policy_year", type = "int", default = 2025L)
tax_increases <- settings$get("figures", "tax_increases", type = "ints", default = 1:5)
states_opt <- settings$get("figures", "states", type = "chrs", default = "all")
make_plots <- settings$get("figures", "plots", type = "bool", default = TRUE)
formats <- settings$get("figures", "formats", type = "chrs", default = "svg")
make_html <- settings$get("figures", "html", type = "bool", default = TRUE) && "svg" %in% formats
workers <- settings$get("figures", "workers", type = "int", default = min(8, detectCores(), na.rm = TRUE))
if (settings$get("grid", "kernel", type = "bool", default = TRUE)) load_kernel(file.path(fast_dir, "lib", "cohort_kernel.cpp"))

prices <- read.csv(prices_csv, stringsAsFactors = FALSE)
prices$State <- toupper(trimws(prices$State))
states <- sort(if (identical(states_opt, "all")) prices$State[prices$State %in% statepop$abbr] else toupper(states_opt),
               method = "radix")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The parts of a runstates() result that the figure scripts read.
keep <- function(out) list(df_mort.outputs = out$df_mort.outputs, df_CSprevs.by.state = out$df_CSprevs.by.state)

# Population weights for the rows of a state's df_CSprevs.by.state: the census population each
# prevalence was averaged over in runstates() (age group x sex x year, in the table's row order).
prev_weights <- function(fipscode) {
  s <- get_state(fipscode, data_dir)
  M <- as.matrix(s$M$census); F <- as.matrix(s$F$census)
  minage <- c(18, 18, 25, 45, 65); maxage <- c(99, 24, 44, 64, 99)
  unlist(lapply(1:5, function(i) {
    rows <- (minage[i] + 1):(maxage[i] + 1)
    dM <- colSums(M[rows, , drop = FALSE]); dF <- colSums(F[rows, , drop = FALSE])
    c(dM, dF, dM + dF)
  }), use.names = FALSE)
}

# The US as all states together: mortality outcomes summed over states, prevalence weighted by
# each state's population. `items` are one scenario's kept results for every state.
us_of <- function(items, wts) {
  mort <- lapply(items, `[[`, "df_mort.outputs")
  prev <- lapply(items, `[[`, "df_CSprevs.by.state")
  for (d in mort[-1]) stopifnot(identical(d[c("year", "gender")], mort[[1]][c("year", "gender")]))
  for (d in prev[-1]) stopifnot(identical(d[c("year", "gender", "age")], prev[[1]][c("year", "gender", "age")]))
  m <- mort[[1]]
  for (cl in setdiff(names(m)[vapply(m, is.numeric, logical(1))], "year")) m[[cl]] <- Reduce(`+`, lapply(mort, `[[`, cl))
  p <- prev[[1]]
  p$prev <- Reduce(`+`, Map(function(d, w) d$prev * w, prev, wts)) / Reduce(`+`, wts)
  m$state <- p$state <- "US"
  m$abbr <- p$abbr <- "US"
  list(df_mort.outputs = m, df_CSprevs.by.state = p)
}

no_effect <- matrix(1, nrow = 100, ncol = cohyears + 100, dimnames = list(NULL, as.character(startbc:endyear)))

t0 <- Sys.time()
res <- mclapply(states, function(st) {
  fipscode <- abbr_to_fips(st)
  initprice <- as.numeric(prices$default_init_price[prices$State == st][1])
  base <- keep(runstates(fipscode, no_effect, no_effect, data_dir = data_dir, full = TRUE))
  taxes <- lapply(tax_increases, function(d) {
    te <- tax_effectCalculation(initprice = initprice, tax = d, startbc = startbc, endyear = endyear,
                                policyYear = policy_year, inidecay = 0.0, cesdecay = 0.2,
                                iniagemod = 1, cesagemod = 1)
    out <- keep(runstates(fipscode, te$m.initiation.effect, te$m.cessation.effect, data_dir = data_dir, full = TRUE))
    c(out, list(initprice = initprice, tax_raise = d))
  })
  names(taxes) <- sprintf("%s_%ddollar", fipscode, tax_increases)
  wts <- prev_weights(fipscode)
  stopifnot(length(wts) == nrow(base$df_CSprevs.by.state))
  list(fips = fipscode, base = base, taxes = taxes, weights = wts)
}, mc.cores = workers)
failed <- vapply(res, function(r) inherits(r, "try-error"), logical(1))
if (any(failed)) stop("model run failed for ", paste(states[failed], collapse = ", "), ": ", as.character(res[[which(failed)[1]]]))

baseline_results <- setNames(lapply(res, `[[`, "base"), vapply(res, `[[`, "", "fips"))
tax_out_list <- unlist(lapply(res, `[[`, "taxes"), recursive = FALSE)
# US entries, keyed "US" and "US_<d>dollar" like the states'.
wts <- lapply(res, `[[`, "weights")
baseline_results <- c(list(US = us_of(lapply(res, `[[`, "base"), wts)), baseline_results)
us_taxes <- lapply(seq_along(tax_increases), function(k) {
  c(us_of(lapply(res, function(r) r$taxes[[k]]), wts), list(initprice = NA_real_, tax_raise = tax_increases[k]))
})
tax_out_list <- c(setNames(us_taxes, sprintf("US_%ddollar", tax_increases)), tax_out_list)
saveRDS(baseline_results, file.path(out_dir, "baseline_results.rds"))
saveRDS(tax_out_list, file.path(out_dir, "tax_out_list.rds"))
message(sprintf("%d states, baseline + %d tax levels from %d: %.1fs -> %s", length(states), length(tax_increases),
                policy_year, as.numeric(difftime(Sys.time(), t0, units = "secs")), out_dir))

if (make_plots) {
  t0 <- Sys.time()
  plots_file <- file.path(repo, "R", "tax_plots.R")
  suppressMessages(source(plots_file))
  fips_vec <- names(baseline_results)
  tax_prevDiff <- tax_prev_diff(fips_vec, baseline_results, tax_out_list, dollar_levels = tax_increases)
  fig_dir <- file.path(out_dir, "figs", "tax")
  # Separate R processes rather than forks: on macOS, graphics devices in a forked child can abort
  # (Objective-C fork safety).
  cl <- makePSOCKcluster(min(workers, length(fips_vec)))
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, c("plots_file", "tax_out_list", "tax_prevDiff", "fig_dir", "tax_increases", "formats"), envir = environment())
  # Every state's prevalence figures also show the US, in grey.
  invisible(clusterEvalQ(cl, suppressMessages(source(plots_file))))
  done <- parLapplyLB(cl, fips_vec, function(s) {
    tryCatch(suppressWarnings(save_tax_plots(s, tax_out_list, tax_prevDiff, fig_dir = fig_dir,
                                             dollar_levels = tax_increases, formats = formats, compare = "US")),
             error = function(e) conditionMessage(e))
  })
  bad <- !vapply(done, is.data.frame, logical(1))
  if (any(bad)) stop("figures failed for ", paste(fips_vec[bad], collapse = ", "), ": ", done[[which(bad)[1]]])
  message(sprintf("figures: %d x %s in %.1fs -> %s", 4 * length(fips_vec), paste(formats, collapse = "+"),
                  as.numeric(difftime(Sys.time(), t0, units = "secs")), fig_dir))
  if (make_html) {
    source(file.path(fast_dir, "lib", "gallery.R"))
    subtitle <- sprintf("Tax increases of %s starting in %d, against no policy. US and %d states, generated %s.",
                        paste0("+$", tax_increases, collapse = ", "), policy_year, length(fips_vec) - 1, format(Sys.Date()))
    page <- write_tax_gallery(file.path(out_dir, "index.html"), do.call(rbind, done), "figs/tax", subtitle)
    message("page: ", page)
  }
}
