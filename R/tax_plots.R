# Tax policy figures: cumulative life-years gained, cumulative deaths avoided, smoking prevalence
# and prevalence reduction by year, for +$1 to +$5 tax increases in one state.
#
# Shared by R/04_visualization_policy_specific.R and R-fast/run_tax_figures.R, so the figures
# are defined in one place. Inputs are the lists 03_main_analysis.R saved:
#   tax_out_list[["<fips>_<d>dollar"]]: runstates() output (df_mort.outputs, df_CSprevs.by.state)
#   baseline_results[["<fips>"]]: runstates() output with no policy
#
# usage: tax_prevDiff <- tax_prev_diff(v_statefips, baseline_results, tax_out_list)
#        save_tax_plots(v_statefips, tax_out_list, tax_prevDiff, fig_dir = "figs/tax")

library(dplyr)
library(ggplot2)
library(scales)

abbr_to_name <- function(ab) {
  names <- c(setNames(state.name, state.abb), DC = "District of Columbia", US = "United States")
  unname(names[ab])
}

# Axis labels: counts as 10k, 1.5M; prevalence ticks on whole percentages.
short_number <- function(x) {
  out <- ifelse(abs(x) >= 1e6, paste0(signif(x / 1e6, 3), "M"),
                ifelse(abs(x) >= 1e3, paste0(signif(x / 1e3, 3), "k"), as.character(signif(x, 3))))
  out[is.na(x)] <- NA
  out
}
whole_breaks <- function(limits) unique(round(scales::breaks_extended(n = 5)(limits)))

# Grey comparison lines (one per tax level), drawn under a figure's coloured lines: used to show
# the US on each state's prevalence figures. NULL (no comparison) adds nothing.
compare_layer <- function(cmp, y) {
  if (is.null(cmp)) return(NULL)
  geom_line(data = cmp, aes(x = year, y = .data[[y]], group = tax_scen), inherit.aes = FALSE,
            colour = "grey72", linewidth = 0.5)
}
compare_caption <- function(caption, cmp, compare) {
  if (is.null(cmp)) caption else sprintf("%s; grey: %s", caption, abbr_to_name(compare))
}

## ----------------------- Tax Visualization -----------------------------------

# ------------------------------------------------------------------------------
# 1. Life Years Gained
# ------------------------------------------------------------------------------

plot_LYG_tax <- function(fipscode,out_list = tax_out_list, dollar_levels = 1:5) {
  df_list <- lapply(dollar_levels, function(d) {
    key <- sprintf("%s_%ddollar", fipscode, d)
    out <- out_list[[key]][["df_mort.outputs"]]
    out$tax_scen <- sprintf("%ddollar", d)
    out
  }) |>bind_rows()
  
  df_plot <- df_list |>
    mutate(gender = recode(tolower(gender), "both" = "Persons")) |>
    filter(gender == "Persons", year >= 2025) |>
    group_by(tax_scen, year) |>
    summarise(LYGcum = sum(LYGcum), .groups = "drop")
  
  state_name <- abbr_to_name(df_list$abbr[1])
  
  tax_cols <- c("1dollar" = "#1b9e77", "2dollar" = "#d95f02", "3dollar" = "#7570b3", "4dollar" = "#e7298a",   "5dollar" = "#66a61e")
        
  ggplot(df_plot, aes(year, LYGcum, colour = tax_scen)) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual("Tax increase", values = tax_cols,
                        labels = paste0("$", dollar_levels)) +
    scale_y_continuous("Life-years Gained", labels = short_number) +
    scale_x_continuous("Year", limits = c(2021, 2100),
                       breaks = seq(2030, 2100, 10)) +
    labs(title = sprintf("%s — Tax Policy", state_name)) +
    theme_light(base_size = 11) +
    theme(legend.position = c(.05, .95),legend.justification = c("left", "top"), legend.box.just = "left",
          legend.background = element_rect(fill = alpha("white", 0.4), colour = NA),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 9),
          axis.text.x = element_text(angle = 60, hjust = 1))}

        
# ------------------------------------------------------------------------------
# 2. Death Avoided
# ------------------------------------------------------------------------------

plot_DeathAvoided_tax <- function(fipscode, out_list = tax_out_list, dollar_levels = 1:5) {
  df_list <- lapply(dollar_levels, function(d) {
    key <- sprintf("%s_%ddollar", fipscode, d)
    out <- out_list[[key]][["df_mort.outputs"]]
    out$tax_scen <- sprintf("%ddollar", d)
    out
  }) |>bind_rows()
  
  df_plot <- df_list |>
    mutate(gender = recode(tolower(gender), "both" = "Persons")) |>
    filter(gender == "Persons", year >= 2025) |>
    group_by(tax_scen, year) |>
    summarise(SADsAvertedcum = sum(SADsAvertedcum), .groups = "drop")
  
  state_name <- abbr_to_name(df_list$abbr[1])
  
  tax_cols <- c("1dollar" = "#1b9e77", "2dollar" = "#d95f02", "3dollar" = "#7570b3", "4dollar" = "#e7298a",   "5dollar" = "#66a61e")
  
  ggplot(df_plot, aes(year, SADsAvertedcum, colour = tax_scen)) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual("Tax increase", values = tax_cols,
                        labels = paste0("$", dollar_levels)) +
    scale_y_continuous("Deaths avoided (cumulative)", labels = short_number) +
    scale_x_continuous("Year", limits = c(2021, 2100),
                       breaks = seq(2030, 2100, 10)) +
    labs(title = sprintf("%s — Tax Policy", state_name)) +
    theme_light(base_size = 11) +
    theme(legend.position = c(.05, .95),legend.justification = c("left", "top"), legend.box.just = "left",
          legend.background = element_rect(fill = alpha("white", 0.4), colour = NA),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 9),
          axis.text.x = element_text(angle = 60, hjust = 1))}


# ------------------------------------------------------------------------------
# 3. Smoking Prevalence
# ------------------------------------------------------------------------------ 

plot_Prev_tax <- function(fipscode, out_list = tax_out_list,dollar_levels = 1:5, compare = NULL) {
  prev_series <- function(code) {
    df_all <- lapply(dollar_levels, function(d) {
      key <- sprintf("%s_%ddollar", code, d)
      df  <- out_list[[key]][["df_CSprevs.by.state"]]
      df$tax_scen <- sprintf("%ddollar", d)
      df}) |> bind_rows()
    df_plot <- df_all |>
      mutate(gender = recode(tolower(gender), "both" = "persons")) |>
      filter(gender == "persons", year >= 2025) |>
      group_by(tax_scen, year) |>
      summarise(prev_pct = mean(prev) * 100, .groups = "drop")
    list(df_all = df_all, df_plot = df_plot)
  }
  series <- prev_series(fipscode)
  df_all <- series$df_all; df_plot <- series$df_plot
  cmp <- if (!is.null(compare) && compare != fipscode) prev_series(compare)$df_plot
  
  state_name <- abbr_to_name(df_all$abbr[1])

  tax_cols <- c("1dollar" = "#1b9e77", "2dollar" = "#d95f02", "3dollar" = "#7570b3", "4dollar" = "#e7298a",   "5dollar" = "#66a61e")
  
  ggplot(df_plot, aes(year, prev_pct, colour = tax_scen)) +
    compare_layer(cmp, "prev_pct") +
    geom_line(linewidth = 0.6) +
    scale_colour_manual("Tax increase", values = tax_cols,
                        labels = paste0("$", dollar_levels)) +
    scale_y_continuous("Smoking prevalence (%)", breaks = whole_breaks, labels = function(x) sprintf("%.0f%%", x)) +
    scale_x_continuous("Year", limits = c(2021, 2100),
                       breaks = seq(2030, 2100, 10)) +
    labs(title = sprintf("%s — Tax Policy", state_name),
         caption = compare_caption("Persons = Both sexes", cmp, compare)) +
    theme_light(base_size = 11) +
    theme(legend.position = c(.85, .95),legend.justification = c("left", "top"), legend.box.just = "left",
          legend.background = element_rect(fill = alpha("white", 0.4), colour = NA),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 9),
          axis.text.x = element_text(angle = 60, hjust = 1)) +
    # With comparison lines, the legend moves above the panel so it does not cover them.
    if (!is.null(cmp)) theme(legend.position = "top", legend.direction = "horizontal",
                             legend.justification = "center", legend.title = element_text(size = 7),
                             legend.text = element_text(size = 7))}

  
# ------------------------------------------------------------------------------
# 4. Prevalence Reduction
# ------------------------------------------------------------------------------

normalise_gender <- function(x){
  tolower(x) |> recode("male" = "Males", "female" = "Females","both" = "Persons")}

# Prevalence reduction (percentage points) of each tax level against the baseline, by state:
# list keyed "<fips>_<d>dollar_diff", as plot_PrevRed_tax() reads it.
tax_prev_diff <- function(fips_vec, baseline_results, tax_out_list, dollar_levels = 1:5) {
  tax_prevDiff <- list()

  for (s in fips_vec){
    base_tbl <- baseline_results[[s]][["df_CSprevs.by.state"]] %>%
      mutate(gender = normalise_gender(gender)) %>%
      select(year, gender, age, prev_base = prev)
  
    for (d in dollar_levels){
      key_tax <- sprintf("%s_%ddollar", s, d)
      pol_tbl <- tax_out_list[[key_tax]][["df_CSprevs.by.state"]] %>%
        mutate(gender = normalise_gender(gender)) %>%
        select(year, gender, age, prev_tax = prev)
    
      diff_tbl <- base_tbl %>%
        inner_join(pol_tbl, by = c("year","gender","age")) %>%
        mutate(prev_diff_pp = (prev_base - prev_tax) * 100) %>%   
        select(year, gender, age, prev_diff_pp)
    
      tax_prevDiff[[paste0(s, "_", d, "dollar_diff")]] <- diff_tbl
    }
  }
  tax_prevDiff
}

plot_PrevRed_tax <- function(fipscode,out_list_diff = tax_prevDiff,  dollar_levels = 1:5, compare = NULL) {
  prevred_series <- function(code) {
    df_all <- lapply(dollar_levels, function(d) {
      key <- sprintf("%s_%ddollar_diff", code, d)
      df  <- out_list_diff[[key]]
      df$tax_scen <- sprintf("%ddollar", d)
      df }) |> bind_rows()
    df_plot <- df_all |>
      mutate(gender = tolower(gender)) |>
      filter(gender == "persons", year >= 2025) |>
      group_by(tax_scen, year) |>
      summarise(prev_red_pp = mean(prev_diff_pp), .groups = "drop")
    list(df_all = df_all, df_plot = df_plot)
  }
  series <- prevred_series(fipscode)
  df_all <- series$df_all; df_plot <- series$df_plot
  cmp <- if (!is.null(compare) && compare != fipscode) prevred_series(compare)$df_plot
 
  
  state_name <- abbr_to_name(df_all$abbr[1])
  
  tax_cols <- c("1dollar" = "#1b9e77", "2dollar" = "#d95f02", "3dollar" = "#7570b3", "4dollar" = "#e7298a",   "5dollar" = "#66a61e")
  
  ggplot(df_plot, aes(year, prev_red_pp, colour = tax_scen)) +
    compare_layer(cmp, "prev_red_pp") +
    geom_line(linewidth = 0.6) +
    scale_colour_manual("Tax increase", values = tax_cols,
                        labels = paste0("$", dollar_levels)) +
    scale_y_continuous("Prevalence reduction (pp)",
                       labels = function(x) sprintf("%.1f pp", x)) +
    scale_x_continuous("Year", limits = c(2021, 2100),
                       breaks = seq(2030, 2100, 10)) +
    labs(title = sprintf("%s — Tax Policy", state_name),
         caption = compare_caption("pp = percentage-points", cmp, compare)) +
    theme_light(base_size = 11) +
    theme(
      legend.position   = "top",             
      legend.direction  = "horizontal",       
      legend.justification = "center",        
      legend.box        = "horizontal",     
      legend.margin     = margin(b = 4),      
      legend.background = element_rect(fill = alpha("white", 0.4), colour = NA),
      legend.title      = element_text(size = 7),
      legend.text       = element_text(size = 7),
      legend.spacing.x  = unit(6, "pt"),  
      legend.key.height = unit(9, "pt"),
      axis.text.x       = element_text(angle = 60, hjust = 1))}

# ------------------------------------------------------------------------------
# Save the four tax figures for each state in fips_vec to fig_dir, as <fips>_<State>_<figure>_tax.<ext>
# for each format in `formats` ("jpg" at 600 dpi, "svg" via svglite, or any ggsave() extension).
# compare: code of entries (e.g. "US") drawn as grey lines on the prevalence figures; NULL for none.
# Returns (invisibly) one row per figure: fips, abbr, figure and the file name of each format.
# ------------------------------------------------------------------------------

save_tax_plots <- function(fips_vec, tax_out_list, tax_prevDiff, fig_dir = "figs/tax", dollar_levels = 1:5,
                           formats = "jpg", compare = NULL) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  plots <- list(LYG = function(s) plot_LYG_tax(s, out_list = tax_out_list, dollar_levels = dollar_levels),
                DeathAvoided = function(s) plot_DeathAvoided_tax(s, out_list = tax_out_list, dollar_levels = dollar_levels),
                Prev = function(s) plot_Prev_tax(s, out_list = tax_out_list, dollar_levels = dollar_levels, compare = compare),
                PrevRed = function(s) plot_PrevRed_tax(s, out_list_diff = tax_prevDiff, dollar_levels = dollar_levels,
                                                   compare = compare))
  written <- list()
  for (s in fips_vec) {
    state_abbr <- tax_out_list[[sprintf("%s_%ddollar", s, dollar_levels[1])]][["df_mort.outputs"]]$abbr[1]
    state_name <- abbr_to_name(state_abbr)
    for (k in names(plots)) {
      p <- plots[[k]](s)
      files <- sprintf("%s_%s_%s_tax.%s", s, gsub(" ", "_", state_name), k, formats)
      for (f in files) ggsave(filename = file.path(fig_dir, f), plot = p, width = 6, height = 4, dpi = 600)
      written[[length(written) + 1]] <- data.frame(fips = s, abbr = state_abbr, figure = k,
                                                   as.list(setNames(files, formats)), stringsAsFactors = FALSE)
    }
  }
  invisible(do.call(rbind, written))
}
