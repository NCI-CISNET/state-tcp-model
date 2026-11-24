# Local test, not needed to produce outputs. Function-level equivalence: runs the original model
# functions (R/02_model_functions.R at a git ref) and lib/tcp_fast.R on the same states and
# scenarios, and requires identical() results for every object they return, including the
# grid runner's tcp_tables() against the original generate_TCPoutput().
#
# usage: Rscript R-fast/tools/check_equivalence.R [--config=<config.yml>]
#          [--states=AK,CA] [--policy-years=2026,2040] [--ref=HEAD]
#   inputs.data_dir and inputs.prices come from the config (or --data-dir / --prices). --ref is the
#   git ref whose R/02_model_functions.R is the reference (default HEAD: this branch's original code).
tools_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])))
fast_dir <- dirname(tools_dir)
repo <- dirname(fast_dir)
source(file.path(fast_dir, "lib", "packages.R"))
use_packages(c("Rcpp", "usmap", "yaml", "dplyr", "readr"))
suppressMessages({ library(dplyr); library(readr) })
source(file.path(fast_dir, "lib", "config.R"))
settings <- read_settings(file.path(fast_dir, "config.yml"))
data_dir <- normalizePath(settings$get("inputs", "data_dir", type = "path", required = TRUE))
prices <- read.csv(settings$get("inputs", "prices", type = "path", required = TRUE), stringsAsFactors = FALSE)
# Test-specific options: command line only, so the grid's states and policy years are not used.
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
states <- strsplit(opt("states", "AK,CA"), ",")[[1]]
policy_years <- as.integer(strsplit(opt("policy-years", "2026,2040"), ",")[[1]])
ref <- opt("ref", "HEAD")

orig_src <- system2("git", c("-C", shQuote(repo), "show", paste0(ref, ":R/02_model_functions.R")), stdout = TRUE)
drop <- grep("^packages <- c\\(|^\\s+\"ggpubr\"|^installed_packages|^if \\(any\\(installed_packages|install\\.packages\\(packages|^lapply\\(packages", orig_src)
orig <- new.env()
eval(parse(text = orig_src[-drop]), envir = orig)
fast <- new.env()
sys.source(file.path(fast_dir, "lib", "tcp_fast.R"), envir = fast)
has_kernel <- fast$load_kernel(file.path(fast_dir, "lib", "cohort_kernel.cpp"))
configs <- list(list(name = "R full", kernel = FALSE, full = TRUE),
                list(name = "R lean", kernel = FALSE, full = FALSE))
if (has_kernel) {
  configs <- c(configs, list(list(name = "kernel full", kernel = TRUE, full = TRUE),
                             list(name = "kernel lean", kernel = TRUE, full = FALSE)))
}

strip_output <- function(l) {
  for (k in c("male_prev_obj", "female_prev_obj", "male_base_prev_obj", "female_base_prev_obj")) l[[k]]$m_output <- NULL
  l
}
check <- function(label, a, b) {
  ok <- identical(a, b)
  cat(sprintf("  %-36s %s\n", label, if (ok) "identical" else "DIFFERENT"))
  if (!ok) print(all.equal(a, b))
  ok
}

prices$State <- toupper(trimws(prices$State))
all_ok <- TRUE
old_wd <- setwd(dirname(data_dir))
for (st in states) {
  fipscode <- fast$abbr_to_fips(st)
  initprice <- as.numeric(prices$default_init_price[prices$State == st][1])
  memo <- new.env()  # tcp_tables() state memo, reused across this state's scenarios as the grid does
  for (tax in c(0.5, 2)) {
    for (py in policy_years) {
      cat(sprintf("%s  base price %.2f  tax +%.2f  policy year %d\n", st, initprice, tax, py))
      te_o <- orig$tax_effectCalculation(initprice = initprice, tax = tax, policyYear = py)
      te_f <- fast$tax_effectCalculation(initprice = initprice, tax = tax, policyYear = py)
      all_ok <- check("tax_effectCalculation", te_o, te_f) && all_ok

      t0 <- Sys.time(); r_o <- orig$runstates(fipscode, te_o$m.initiation.effect, te_o$m.cessation.effect)
      t_o <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      g_o <- orig$generate_TCPoutput(r_o$l_prev_out, fipscode, py, policy_name = "taxes", write_files = FALSE,
                                     out_root = file.path(tempdir(), "tcp_tool"))
      times <- c(original = t_o)

      for (cfg in configs) {
        options(tcp_fast.kernel = cfg$kernel)
        fast$.cache$key <- NULL
        invisible(fast$get_state(fipscode, data_dir))
        t0 <- Sys.time()
        r_f <- fast$runstates(fipscode, te_f$m.initiation.effect, te_f$m.cessation.effect,
                              data_dir = data_dir, full = cfg$full)
        times[cfg$name] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        if (cfg$full) {
          all_ok <- check(paste(cfg$name, "l_prev_out"), strip_output(r_o$l_prev_out), strip_output(r_f$l_prev_out)) && all_ok
          all_ok <- check(paste(cfg$name, "df_CSprevs.by.state"), r_o$df_CSprevs.by.state, r_f$df_CSprevs.by.state) && all_ok
          all_ok <- check(paste(cfg$name, "df_mort.outputs"), r_o$df_mort.outputs, r_f$df_mort.outputs) && all_ok
        }
        g_f <- fast$generate_TCPoutput(r_f$l_prev_out, fipscode, py, policy_name = "taxes", write_files = FALSE)
        for (k in c("deaths", "lyg", "results")) {
          a <- g_o[[k]]; b <- g_f[[k]]
          rownames(a) <- rownames(b) <- NULL
          all_ok <- check(paste(cfg$name, "TCP", k), a, b) && all_ok
        }
        if (!cfg$full) {
          tt <- fast$tcp_tables(r_f$l_prev_out, py, memo)
          for (k in c("deaths", "lyg", "results")) {
            a <- g_o[[k]]; rownames(a) <- NULL
            all_ok <- check(paste(cfg$name, "tcp_tables", k), a, fast$tcp_table_df(tt[[k]])) && all_ok
          }
        }
      }
      cat("  runstates time:", paste(sprintf("%s %.3fs", names(times), times), collapse = ", "), "\n")
    }
  }
}
setwd(old_wd)
cat(if (all_ok) "\nALL CHECKS IDENTICAL\n" else "\nSOME CHECKS DIFFER\n")
quit(status = if (all_ok) 0 else 1)
