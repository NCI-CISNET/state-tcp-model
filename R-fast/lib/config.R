# Settings for the R-fast scripts: a YAML config file plus command-line overrides.
#
# Every script reads R-fast/config.yml, or the file given with --config=<path>. A
# command-line flag --<name>=<value> overrides the setting of that name, so a one-off run does
# not need an edited config. Relative paths in the config file are relative to the config file
# itself (so it can live anywhere); relative paths given on the command line are relative to the
# working directory. config.example.yml lists every setting.
#
#   s <- read_settings(file.path(fast_dir, "config.yml"))
#   data_dir <- s$get("inputs", "data_dir", type = "path", required = TRUE)
#   workers  <- s$get("grid", "workers", type = "int", default = 8)   # --workers=N overrides

read_settings <- function(default_config, args = commandArgs(trailingOnly = TRUE)) {
  flag_value <- function(name) {
    hit <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else NULL
  }
  given <- flag_value("config")
  cfg_path <- if (is.null(given)) default_config else given
  cfg <- list()
  cfg_dir <- getwd()
  if (file.exists(cfg_path)) {
    cfg_path <- normalizePath(cfg_path)
    cfg <- yaml::read_yaml(cfg_path)
    cfg_dir <- dirname(cfg_path)
  } else if (!is.null(given)) {
    stop("Config file not found: ", cfg_path)
  } else {
    cfg_path <- NA_character_
  }

  # section/key: where the setting sits in the YAML file. flag: its command-line name (key with
  # "_" as "-" unless given). type: path, chr, chrs, int, ints, nums or bool.
  get <- function(section, key, type = "chr", default = NULL, flag = gsub("_", "-", key), required = FALSE) {
    v <- flag_value(flag)
    from_cli <- !is.null(v)
    if (!from_cli) v <- cfg[[section]][[key]]
    if (is.null(v) || (length(v) == 1 && is.character(v) && !nzchar(v))) {
      if (required) {
        stop(sprintf("Setting %s.%s is required: set it in %s or pass --%s=<value>",
                     section, key, if (is.na(cfg_path)) "a config file (see config.example.yml)" else cfg_path, flag))
      }
      return(default)
    }
    split <- function(x) if (from_cli) trimws(strsplit(x, ",")[[1]]) else unlist(x)
    switch(type,
      path = {
        p <- path.expand(as.character(v))
        if (!from_cli && !grepl("^/", p)) p <- file.path(cfg_dir, p)
        p
      },
      # YAML reads on/off/yes/no as logicals; settings such as link = off want the word.
      chr = if (is.logical(v)) (if (v) "on" else "off") else as.character(v),
      chrs = as.character(split(v)),
      int = as.integer(v),
      ints = as.integer(split(v)),
      nums = as.numeric(split(v)),
      bool = if (is.logical(v)) v else tolower(v) %in% c("on", "true", "yes", "1"),
      stop("unknown setting type: ", type))
  }
  list(get = get, file = cfg_path)
}
