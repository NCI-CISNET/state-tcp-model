# Local test, not needed to produce outputs. Compares two TCP tool output trees
# (<ST>/<policy year>/<type>/<file>.csv): byte-identical files, and for files that differ, how many
# cells differ and by how much. Symlinked files are compared by their contents.
#
# usage: Rscript R-fast/tools/compare_trees.R --ref=<dir> --new=<dir> [--states=AK,CA,...] [--workers=N]
suppressMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
ref <- normalizePath(opt("ref")); new <- normalizePath(opt("new"))
workers <- as.integer(opt("workers", detectCores()))
tree_files <- function(root) {
  f <- list.files(root, pattern = "\\.csv$", recursive = TRUE)
  f <- f[lengths(strsplit(f, "/")) == 4]
  states <- opt("states")
  if (!is.null(states)) f <- f[sub("/.*", "", f) %in% c(strsplit(states, ",")[[1]], "US")]
  f
}
fa <- tree_files(ref); fb <- tree_files(new)
only <- c(setdiff(fa, fb), setdiff(fb, fa))
common <- intersect(fa, fb)

compare_file <- function(f) {
  parts <- strsplit(f, "/")[[1]]
  kind <- paste(if (parts[1] == "US") "US" else "state", parts[3])
  a <- readBin(file.path(ref, f), "raw", file.size(file.path(ref, f)))
  b <- readBin(file.path(new, f), "raw", file.size(file.path(new, f)))
  if (identical(a, b)) return(data.frame(kind = kind, status = "identical", cells = 0, of = 0, abs = 0, rel = 0))
  x <- read.csv(file.path(ref, f), colClasses = "character", na.strings = character())
  y <- read.csv(file.path(new, f), colClasses = "character", na.strings = character())
  if (!identical(dim(x), dim(y)) || !identical(names(x), names(y))) {
    return(data.frame(kind = kind, status = "shape differs", cells = NA, of = NA, abs = Inf, rel = Inf))
  }
  nd <- 0; worst <- 0; worst_rel <- 0
  for (cl in names(x)) {
    m <- x[[cl]] != y[[cl]]
    if (!any(m)) next
    nd <- nd + sum(m)
    xa <- as.numeric(x[[cl]][m]); ya <- as.numeric(y[[cl]][m])
    d <- abs(xa - ya)
    worst <- max(worst, d, na.rm = TRUE)
    worst_rel <- max(worst_rel, d / pmax(abs(xa), abs(ya)), na.rm = TRUE)
  }
  data.frame(kind = kind, status = "differs", cells = nd, of = prod(dim(x)), abs = worst, rel = worst_rel)
}

rows <- do.call(rbind, mclapply(common, compare_file, mc.cores = workers, mc.preschedule = TRUE))
summary <- do.call(rbind, lapply(split(rows, list(rows$kind, rows$status), drop = TRUE), function(d) {
  data.frame(type = d$kind[1], status = d$status[1], files = nrow(d), cells_differ = sum(d$cells),
             of_cells = sum(d$of), max_abs_diff = max(d$abs), max_rel_diff = max(d$rel))
}))
rownames(summary) <- NULL
cat(sprintf("=== %s vs %s: %d files ===\n", new, ref, length(fa)))
print(summary[order(summary$type, summary$status), ], row.names = FALSE)
if (length(only)) cat(sprintf("files present in only one tree: %d (e.g. %s)\n", length(only), paste(head(only, 3), collapse = ", ")))
ok <- length(only) == 0 && all(rows$status == "identical")
cat(if (ok) "\nALL IDENTICAL\n" else "\nDIFFERENCES FOUND\n")
quit(status = if (ok) 0 else 1)
