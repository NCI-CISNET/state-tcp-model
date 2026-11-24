# TCP tool files for the cigarette-tax scenario grid, for every state and the US, computed in
# parallel with the fast model in lib/tcp_fast.R. See README.md for what is computed and how it
# differs from the original R/ pipeline.
#
# The grid is every price-per-pack target (grid.ppp_targets) crossed with every added tax
# (grid.added_taxes), for each policy year (grid.policy_years). A scenario's effect depends only on
# its total tax increase, max(0, target - state price) + added tax, so tags with the same total
# give identical outputs: only one representative per total is run, and jobs.csv records which
# file tags each representative covers. Tags that share a representative link to one file rather
# than holding copies: link = symlink (the default) makes relative symlinks, the cheapest on APFS;
# hard makes hard links; off writes copies. aws s3 cp/sync follows symlinks, so a bucket gets
# real files either way.
#
# grid.final_dir: the TCP tool files that processing.ipynb and update_us_results.py produced from
#   the original per-scenario outputs: <ST>/<policy year>/{deaths,lyg,results}/<type>_<ppp>_t<tax>.csv
#   for every tag, with *_diff columns in results and values at 6 significant digits; and the same
#   for US: deaths and lyg summed over the states run, and (with inputs.population_csv) results
#   weighted by the standardized model population. The US files are built from what the workers
#   return, without reading any files back. Values are the model's exact doubles; the notebook's
#   were read back from CSV by pandas, which left rounding noise in a few near-zero US sums (on the
#   2025 inputs, 7 of 98,280 files have one cell that is ~1e-29 there and exactly 0 here).
# grid.raw_dir: one CSV per (state, representative scenario, output type) with all policy years
#   stacked, as the original per-scenario runs wrote them: <ST>_{deaths,lyg,results}_<ppp>_t<tax>.csv
#
# Scenario cache (grid.cache, default on; final files only): each representative's results are
# saved under a key made from everything they depend on: the state's four input files (md5), the
# initial price and total tax, the policy years, the model code and the population file. A later
# run with the same key reuses the saved results instead of running the model; any change to those
# gives a new key, so the scenario is recomputed. Cache files live in
# tools::R_user_dir("tcp_fast", "cache")/scenarios and are not removed automatically.
#
# usage: Rscript R-fast/run_tax_grid.R [--config=<config.yml>] [--<setting>=<value> ...]
#   Settings come from R-fast/config.yml (see config.example.yml); a flag overrides one, e.g.
#   --states=AK,CA --final-dir=/tmp/out --cache=off. Settings: data-dir, prices, population-csv,
#   final-dir, raw-dir, states, policy-years, ppp-targets, added-taxes, workers, chunk-size,
#   kernel, link, cache.
fast_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])))
source(file.path(fast_dir, "lib", "packages.R"))
use_packages(c("Rcpp", "readr", "vroom", "usmap", "yaml"))
suppressMessages({ library(readr); library(parallel) })
source(file.path(fast_dir, "lib", "config.R"))
source(file.path(fast_dir, "lib", "tcp_fast.R"))
source(file.path(fast_dir, "lib", "us_results.R"))

settings <- read_settings(file.path(fast_dir, "config.yml"))
data_dir <- normalizePath(settings$get("inputs", "data_dir", type = "path", required = TRUE))
prices_csv <- settings$get("inputs", "prices", type = "path", required = TRUE)
population_csv <- settings$get("inputs", "population_csv", type = "path")
raw_dir <- settings$get("grid", "raw_dir", type = "path")
final_dir <- settings$get("grid", "final_dir", type = "path")
if (is.null(raw_dir) && is.null(final_dir)) stop("Set grid.final_dir and/or grid.raw_dir (or --final-dir / --raw-dir)")
policy_years <- settings$get("grid", "policy_years", type = "ints", default = c(2026:2030, 2035, 2040))
ppp_targets <- settings$get("grid", "ppp_targets", type = "nums", default = seq(8, 15, by = 0.5))
added_taxes <- settings$get("grid", "added_taxes", type = "nums", default = c(0, 1, 1.5, 2, 2.5, 3))
states_opt <- settings$get("grid", "states", type = "chrs", default = "all")
# Throughput on a 6 performance + 6 efficiency core M3 Pro peaks at about 8 workers; more adds
# memory contention without adding speed.
workers <- settings$get("grid", "workers", type = "int", default = min(8, detectCores(), na.rm = TRUE))
chunk_size <- settings$get("grid", "chunk_size", type = "int", default = 0L)
link_mode <- settings$get("grid", "link", default = "symlink")
if (link_mode == "on") link_mode <- "hard"  # the old --link=on
if (!link_mode %in% c("symlink", "hard", "off")) stop("link must be symlink, hard or off")
# Raw mode always computes: the raw files are the reference the final files are checked against.
use_cache <- settings$get("grid", "cache", type = "bool", default = TRUE) && is.null(raw_dir) && !is.null(final_dir)
if (settings$get("grid", "kernel", type = "bool", default = TRUE)) load_kernel(file.path(fast_dir, "lib", "cohort_kernel.cpp"))
value_cols <- list(deaths = c("deaths_avoided_males", "deaths_avoided_females"),
                   lyg = c("cLYG_males", "cLYG_females"),
                   results = c("males_baseline", "females_baseline", "males_policy", "females_policy",
                               "both_baseline", "both_policy"))

prices <- read.csv(prices_csv, stringsAsFactors = FALSE)
prices$State <- toupper(trimws(prices$State))
states <- sort(if (identical(states_opt, "all")) prices$State[prices$State %in% statepop$abbr] else toupper(states_opt),
               method = "radix")

build_jobs <- function(st) {
  initprice <- as.numeric(prices$default_init_price[prices$State == st][1])
  jobs <- list()
  for (tax in added_taxes) {
    for (ppp in ppp_targets) {
      key <- sprintf("%.4f", max(0, ppp - initprice) + tax)
      tag <- sprintf("%0.2f_t%0.2f", ppp, tax)
      if (is.null(jobs[[key]])) jobs[[key]] <- list(state = st, initprice = initprice, ppp = ppp, tax = tax, tags = tag)
      else jobs[[key]]$tags <- c(jobs[[key]]$tags, tag)
    }
  }
  unname(jobs)
}
jobs <- unlist(lapply(states, build_jobs), recursive = FALSE)
rep_tag <- function(job) sprintf("%0.2f_t%0.2f", job$ppp, job$tax)

jobs_df <- data.frame(state = sapply(jobs, `[[`, "state"), ppp = sapply(jobs, `[[`, "ppp"),
                      tax = sapply(jobs, `[[`, "tax"), tags = sapply(jobs, function(j) paste(j$tags, collapse = ";")))
for (d in c(raw_dir, final_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  write.csv(jobs_df, file.path(d, "jobs.csv"), row.names = FALSE)
}
if (!is.null(final_dir)) {
  for (st in c(states, "US")) for (py in policy_years) for (k in names(value_cols)) {
    dir.create(file.path(final_dir, st, py, k), recursive = TRUE, showWarnings = FALSE)
  }
}

# Hash of a string (written to a file, then md5), and of a set of files.
text_hash <- function(s) {
  f <- tempfile()
  on.exit(unlink(f))
  writeLines(s, f)
  unname(tools::md5sum(f))
}

# Scenario cache keys: everything a representative's results depend on.
cache_root <- file.path(tools::R_user_dir("tcp_fast", "cache"), "scenarios")
code_hash <- pop_hash <- state_hash <- NULL
if (use_cache) {
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  code_files <- c(file.path(fast_dir, "run_tax_grid.R"),
                  file.path(fast_dir, "lib", c("tcp_fast.R", "us_results.R", "cohort_kernel.cpp")))
  code_hash <- text_hash(paste(unname(tools::md5sum(code_files)), collapse = " "))
  pop_hash <- if (is.null(population_csv)) "none" else text_hash(unname(tools::md5sum(population_csv)))
  state_files <- function(st) {
    f <- abbr_to_fips(st)
    file.path(data_dir, "state_inputs", c(file.path("mort_rates", paste0("p.mort_", f, ".RData")),
                                          paste0("smk_", f, ".RData"), paste0("pop_", f, ".RData"),
                                          paste0("le_", f, ".RData")))
  }
  t0 <- Sys.time()
  state_hash <- setNames(unlist(mclapply(states, function(st) {
    text_hash(paste(unname(tools::md5sum(state_files(st))), collapse = " "))
  }, mc.cores = workers)), states)
  message(sprintf("input hashes: %.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
scenario_key <- function(job) {
  total_tax <- max(0, job$ppp - job$initprice) + job$tax
  text_hash(paste(job$state, state_hash[[job$state]], sprintf("%.10f", job$initprice), sprintf("%.10f", total_tax),
                  paste(policy_years, collapse = ","), code_hash, pop_hash, sep = "|"))
}

# Loaded before the workers fork so they share it.
us_pop <- NULL
if (!is.null(final_dir) && !is.null(population_csv)) {
  t0 <- Sys.time()
  us_pop <- load_us_population(population_csv)
  message(sprintf("population weights: %.1fs", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

# Tables are list(labels = <data frame of the non-value columns, including policy_year>,
# values = <matrix of the value columns>), stacked over policy years (see tcp_tables()).

fmt6 <- function(x) {
  s <- sprintf("%.6g", x)
  s[is.na(x)] <- ""
  s
}

# One policy year of one output type, formatted as processing.ipynb wrote it with pandas. The
# label text and, for results, the baseline columns are the same for every scenario of a state, so
# they are formatted once and kept in `fmt` (an environment per state).
format_final <- function(tab, type, py, fmt) {
  key <- paste(type, py)
  f <- fmt[[key]]
  if (is.null(f)) {
    rows <- which(tab$labels$policy_year == py)
    l <- tab$labels[rows, , drop = FALSE]
    vc <- colnames(tab$values)
    if (type == "results") {
      f <- list(rows = rows, pre = paste(l$year, l$age, l$cohort, sep = ","),
                header = paste(c("year", "age", "cohort", vc, "policy_year", "males_diff", "females_diff", "both_diff"),
                               collapse = ","))
    } else {
      f <- list(rows = rows, pre = paste(l$year, l$cohort, sep = ","),
                header = paste(c("year", "cohort", vc, "policy_year"), collapse = ","))
    }
    f$post <- as.character(py)
    assign(key, f, envir = fmt)
  }
  v <- tab$values[f$rows, , drop = FALSE]
  if (type == "results") {
    if (is.null(f$base)) {  # baseline columns: the same for every scenario (and US tag)
      f$base <- lapply(c(1, 2, 5), function(j) fmt6(v[, j]))
      assign(key, f, envir = fmt)
    }
    base <- f$base
    lines <- paste(f$pre, base[[1]], base[[2]], fmt6(v[, 3]), fmt6(v[, 4]), base[[3]], fmt6(v[, 6]), f$post,
                   fmt6(v[, 1] - v[, 3]), fmt6(v[, 2] - v[, 4]), fmt6(v[, 5] - v[, 6]), sep = ",")
  } else {
    lines <- paste(f$pre, fmt6(v[, 1]), fmt6(v[, 2]), f$post, sep = ",")
  }
  c(f$header, lines)
}

# Tags with the same content share one file: the first is written, the others link to it
# (relative symlinks in the same directory, or hard links).
write_final <- function(tables, st, tags, fmt) {
  for (k in names(tables)) {
    for (py in policy_years) {
      txt <- format_final(tables[[k]], k, py, fmt)
      paths <- file.path(final_dir, st, py, k, paste0(k, "_", tags, ".csv"))
      # A link shares its file with the others, so writing in place would change them all.
      # expand = FALSE: the paths are literal, and wildcard expansion costs a directory scan each.
      unlink(paths, expand = FALSE)
      writeLines(txt, paths[1], useBytes = TRUE)
      for (p in paths[-1]) {
        linked <- switch(link_mode, symlink = file.symlink(basename(paths[1]), p),
                         hard = file.link(paths[1], p), off = FALSE)
        if (!linked) writeLines(txt, p, useBytes = TRUE)
      }
    }
  }
}

# Runs the model for one representative scenario. Returns the raw tables (with --raw-dir), the final
# tables (with --final-dir) and what the US files need from it. `cache` is the state's
# environment: its memo holds what tcp_tables() computes once per state and policy year.
compute_scenario <- function(job, cache) {
  fipscode <- abbr_to_fips(job$state)
  total_tax <- max(0, job$ppp - job$initprice) + job$tax
  if (is.null(cache$memo)) cache$memo <- new.env()
  per_py <- lapply(policy_years, function(py) {
    te <- tax_effectCalculation(initprice = job$initprice, tax = total_tax,
                                startbc = startbc, endyear = endyear, policyYear = py,
                                inidecay = 0.0, cesdecay = 0.2, iniagemod = 1, cesagemod = 1)
    out <- runstates(fipscode, te$m.initiation.effect, te$m.cessation.effect, data_dir = data_dir,
                     full = FALSE)
    tcp_tables(out$l_prev_out, py, cache$memo)
  })
  types <- setNames(nm = names(value_cols))
  if (is.null(cache$labels)) {
    cache$labels <- lapply(types, function(k) {
      d <- do.call(rbind, lapply(per_py, function(t) t[[k]]$labels))
      rownames(d) <- NULL
      d
    })
  }
  tables <- lapply(types, function(k) {
    list(labels = cache$labels[[k]], values = do.call(rbind, lapply(per_py, function(t) t[[k]]$values)))
  })
  raw <- if (!is.null(raw_dir)) lapply(tables, tcp_table_df) else NULL
  if (is.null(final_dir)) return(list(raw = raw, final = NULL, part = NULL))
  part <- list(state = job$state, tags = job$tags, deaths = tables$deaths$values, lyg = tables$lyg$values)
  if (!is.null(us_pop)) {
    if (is.null(cache$weights)) cache$weights <- us_row_weights(us_pop[[job$state]], tables$results$labels)
    # The values update_us_results.py reads back from this state's 6-significant-digit files.
    v <- tables$results$values
    vals <- lapply(setNames(nm = value_cols$results[1:4]), function(cl) signif(v[, cl], 6))
    part$us <- us_terms(vals, cache$weights)
  }
  list(raw = raw, final = tables, part = part)
}

# Writes one representative's files, from the cache when its key matches. Returns its US inputs
# (part) and whether it was a cache hit.
run_scenario <- function(job, cache) {
  cache_file <- NULL
  saved <- NULL
  if (use_cache) {
    cache_file <- file.path(cache_root, paste0(scenario_key(job), ".rds"))
    if (file.exists(cache_file)) saved <- readRDS(cache_file)
  }
  hit <- !is.null(saved)
  if (!hit) {
    saved <- compute_scenario(job, cache)
    if (!is.null(cache_file)) {
      # Written under a temporary name and renamed, so a half-written file is never read.
      tmp <- paste0(cache_file, ".", Sys.getpid(), ".tmp")
      saveRDS(saved, tmp)
      file.rename(tmp, cache_file)
    }
  }
  if (!is.null(raw_dir)) {
    for (k in names(saved$raw)) {
      write_csv(saved$raw[[k]], file.path(raw_dir, paste0(job$state, "_", k, "_", rep_tag(job), ".csv")), num_threads = 1)
    }
  }
  if (!is.null(final_dir)) {
    if (is.null(cache$fmt)) cache$fmt <- new.env()
    write_final(saved$final, job$state, job$tags, cache$fmt)
  }
  list(part = saved$part, hit = hit,
       labels = if (!is.null(saved$final)) lapply(saved$final, `[[`, "labels"))
}

# Each chunk holds scenarios from a single state so the state's inputs, baseline and population
# weights are computed once per chunk. By default chunks are sized to give each worker about three,
# and the largest are scheduled first to balance the workers.
if (chunk_size <= 0) chunk_size <- max(1, ceiling(length(jobs) / (workers * 3)))
by_state <- split(jobs, sapply(jobs, `[[`, "state"))
chunks <- unlist(lapply(by_state, function(js) {
  split(js, ceiling(seq_along(js) / chunk_size))
}), recursive = FALSE)
chunks <- chunks[order(-lengths(chunks))]

message(sprintf("%d states, %d scenarios x %d policy years, %d chunks, %d workers, %s, cache %s, links %s",
                length(states), length(jobs), length(policy_years), length(chunks), workers,
                if (use_kernel()) "compiled kernel" else "R implementation",
                if (use_cache) "on" else "off", link_mode))
t_start <- Sys.time()
res <- mclapply(chunks, function(chunk) {
  t0 <- Sys.time()
  cache <- new.env()
  parts <- list()
  layout <- NULL
  hits <- 0
  for (job in chunk) {
    r <- run_scenario(job, cache)
    hits <- hits + r$hit
    p <- r$part
    if (is.null(p)) next
    if (is.null(layout)) layout <- r$labels
    else if (!identical(layout, r$labels)) stop("output layout differs within ", job$state)
    parts[[length(parts) + 1]] <- p
  }
  list(n = length(chunk), secs = as.numeric(difftime(Sys.time(), t0, units = "secs")),
       parts = parts, layout = layout, hits = hits)
}, mc.cores = workers, mc.preschedule = FALSE)

failed <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
if (any(failed)) {
  for (i in which(failed)) message("FAILED chunk ", names(chunks)[i], ": ", as.character(res[[i]]))
  quit(status = 1)
}
wall <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
busy <- sum(sapply(res, `[[`, "secs"))
message(sprintf("model + state files: %d scenarios in %.1fs wall, %.1f worker-seconds (%.2fs per scenario per worker), %.1f scenarios/s, %d from cache",
                length(jobs), wall, busy, busy / length(jobs), length(jobs) / wall, sum(sapply(res, `[[`, "hits"))))

if (!is.null(final_dir)) {
  # US files, from the workers' returned values: deaths and lyg are the state tables added in
  # sorted state order with missing values counting as 0 (DataFrame.add(fill_value = 0) in
  # processing.ipynb); results are population weighted as update_us_results.py does.
  t0 <- Sys.time()
  layout <- res[[1]]$layout
  for (r in res) if (!identical(r$layout, layout)) stop("output layout differs between states")
  parts <- unlist(lapply(res, `[[`, "parts"), recursive = FALSE)
  rm(res)
  part_of <- new.env()
  for (i in seq_along(parts)) for (tag in parts[[i]]$tags) assign(paste(parts[[i]]$state, tag), i, envir = part_of)
  tags <- sort(unique(unlist(lapply(parts, `[[`, "tags"))), method = "radix")
  add_fill <- function(a, b) {
    s <- a + b
    s[is.na(a) & !is.na(b)] <- b[is.na(a) & !is.na(b)]
    s[!is.na(a) & is.na(b)] <- a[!is.na(a) & is.na(b)]
    s
  }
  us_fmt <- new.env()  # label and baseline text, shared by every US tag (each worker keeps a copy)
  us_done <- mclapply(tags, function(tag) {
    idx <- vapply(states, function(st) get(paste(st, tag), envir = part_of), integer(1))
    out <- list()
    for (k in c("deaths", "lyg")) {
      m <- parts[[idx[1]]][[k]]
      for (i in idx[-1]) m <- add_fill(m, parts[[i]][[k]])
      out[[k]] <- list(labels = layout[[k]], values = m)
    }
    if (!is.null(us_pop)) {
      us <- us_results_table(us_fold(lapply(idx, function(i) parts[[i]]$us)), layout$results)
      out$results <- list(labels = us[c("year", "age", "cohort", "policy_year")],
                          values = as.matrix(us[value_cols$results]))
    }
    write_final(out, "US", tag, us_fmt)
    TRUE
  }, mc.cores = workers)
  bad <- !vapply(us_done, isTRUE, logical(1))
  if (any(bad)) stop("US files failed: ", as.character(us_done[which(bad)[1]]))
  message(sprintf("US files (%s): %.1fs", if (is.null(us_pop)) "deaths, lyg; no --population-csv so no results" else "deaths, lyg, results",
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
