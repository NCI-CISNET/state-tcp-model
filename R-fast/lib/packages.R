# Makes sure the R packages a script needs are installed: any that are missing are installed from
# CRAN the first time the script runs. Packages are not attached; scripts call library() for the
# ones they use unqualified.
#
#   use_packages(c("Rcpp", "readr", "yaml"))

use_packages <- function(pkgs) {
  installed <- function(p) vapply(p, requireNamespace, logical(1), quietly = TRUE)
  missing <- pkgs[!installed(pkgs)]
  if (!length(missing)) return(invisible(TRUE))
  repos <- getOption("repos")
  if (is.null(repos) || !nzchar(repos["CRAN"]) || identical(unname(repos["CRAN"]), "@CRAN@")) {
    repos <- c(CRAN = "https://cloud.r-project.org")
  }
  message("Installing missing R packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = repos)
  still <- missing[!installed(missing)]
  if (length(still)) stop("Could not install R packages: ", paste(still, collapse = ", "))
  invisible(TRUE)
}
