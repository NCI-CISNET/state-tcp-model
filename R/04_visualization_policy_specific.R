
# Policy Effects Visualization: LYG, Death Avoided, Smoking Prevalence, Prevalence Reduction

# Plot Function for tax: plot_LYG_tax(),plot_DeathAvoided_tax(), plot_Prev_tax(), plot_PrevRed_tax()
# Plot Function for airlaws: plot_LYG_airlaws(),plot_DeathAvoided_airlaws(), plot_Prev_airlaws(), plot_PrevRed_airlaws()
# Plot Function for tcexp: plot_LYG_tcexp(),plot_DeathAvoided_tcexp(), plot_Prev_tcexp(), plot_PrevRed_tcexp()
# Plot Function for warnings: plot_LYG_warnings(),plot_DeathAvoided_warnings(), plot_Prev_warnings(), plot_PrevRed_warnings()

# PrevRed list: tax_prevDiff(),airlaws_prevDiff(), tcexp_prevDiff(), warnings_prevDiff()

library(dplyr)
library(ggplot2)
library(scales)
library(fs)
library(purrr)


## ----------------------- Tax Visualization -----------------------------------

# The tax figure functions (and abbr_to_name) live in R/tax_plots.R, shared with
# R-fast/run_tax_figures.R.
source("R/tax_plots.R")

tax_prevDiff <- tax_prev_diff(v_statefips, baseline_results, tax_out_list)
save_tax_plots(v_statefips, tax_out_list, tax_prevDiff, fig_dir = "figs/tax")

## ----------------------- Smoke-free Air Laws Visualization -------------------

# ------------------------------------------------------------------------------
# 1. Life Years Gained
# ------------------------------------------------------------------------------
plot_LYG_airlaws <- function(fipscode, scenarios = NULL, out_list = airlaws_out_list) {
  if (is.null(scenarios)) {
    scenarios <- apply(expand.grid(Iwp = 0:1, Ir = 0:1, Ib = 0:1), 1,
                       \(v) sprintf("wp%d_r%d_b%d", v[1], v[2], v[3]))
    scenarios <- scenarios[scenarios != "wp0_r0_b0"]}
  
  df_all <- map_dfr(scenarios, \(scen) {
    key <- sprintf("%s_air_%s", fipscode, scen)
    df <- out_list[[key]][["df_mort.outputs"]]
    df$scenario <- scen
    df})
  
  df_plot <- df_all %>%
    filter(tolower(gender) == "both", year >= 2025) %>%
    group_by(year, scenario) %>%
    summarise(LYGcum = sum(LYGcum), .groups = "drop") %>%
    group_by(scenario) %>%
    filter(any(LYGcum != 0) | scenario == "wp1_r1_b1") %>%
    ungroup()
  
  levels_ordered <- c(
    "wp0_r0_b1", "wp0_r1_b0", "wp0_r1_b1",
    "wp1_r0_b0", "wp1_r0_b1", "wp1_r1_b0", "wp1_r1_b1")
  
  df_plot <- df_plot %>%
    filter(scenario %in% levels_ordered) %>%
    mutate(scenario = factor(scenario, levels = levels_ordered))
  
  state_name <- abbr_to_name(df_all$abbr[1])
  
  ggplot(df_plot, aes(year, LYGcum, colour = scenario)) +
    geom_line(linewidth = 0.6) + 
    scale_y_continuous("Life-years gained", limits = c(0, NA), labels = comma) +
    scale_x_continuous("Year", limits = c(2021, 2100), breaks = seq(2030, 2100, 10)) +
    labs(title  = sprintf("%s — Smoke-free Air-Law Scenarios", state_name), colour = "Scenario") +
    theme_light(base_size = 11) + 
    theme(legend.position = c(.05, .95), legend.justification = c("left", "top"), legend.box.just = "left",
          legend.background = element_rect(fill = alpha("white", .4), colour = NA),
          legend.title = element_text(size = 9), legend.text = element_text(size = 9), axis.text.x = element_text(angle = 60, hjust = 1))}

for (s in v_statefips) {
  sample_key <- grep(paste0("^", s, "_air_"), names(airlaws_out_list), value = TRUE)[1]
  abbr <- airlaws_out_list[[sample_key]][["df_mort.outputs"]]$abbr[1]
  state_name <- abbr_to_name(abbr)
  
  ggsave(filename = sprintf("figs/airlaws/%s_%s_LYG_airlaws.jpg", s, gsub(" ", "_", state_name)),
         plot = plot_LYG_airlaws(s),
         width = 6, height = 4, dpi = 600)}

# ------------------------------------------------------------------------------
# 2. Death Avoided
# ------------------------------------------------------------------------------
plot_DeathAvoided_airlaws <- function(fipscode, scenarios = NULL, out_list = airlaws_out_list) {
  if (is.null(scenarios)) {
    scenarios <- apply(expand.grid(Iwp = 0:1, Ir = 0:1, Ib = 0:1), 1,
                       \(v) sprintf("wp%d_r%d_b%d", v[1], v[2], v[3]))
    scenarios <- scenarios[scenarios != "wp0_r0_b0"]}
  
  df_all <- lapply(scenarios, function(scen) {
    key <- sprintf("%s_air_%s", fipscode, scen)
    df <- out_list[[key]][["df_mort.outputs"]]
    df$scenario <- scen
    df }) |> bind_rows()
  
  df_plot <- df_all |>
    filter(tolower(gender) == "both", year >= 2025) |>
    group_by(year, scenario) |>
    summarise(SADsAvertedcum = sum(SADsAvertedcum), .groups = "drop") |>
    group_by(scenario) |>
    filter(any(SADsAvertedcum != 0) | scenario == "wp1_r1_b1") |>
    ungroup()
  
  scenario_levels <- c(
    "wp0_r0_b1", "wp0_r1_b0", "wp0_r1_b1",
    "wp1_r0_b0", "wp1_r0_b1", "wp1_r1_b0", "wp1_r1_b1")
  
  df_plot <- df_plot |>
    filter(scenario %in% scenario_levels) |>
    mutate(scenario = factor(scenario, levels = scenario_levels))
  
  state_name <- abbr_to_name(df_all$abbr[1])
  
  ggplot(df_plot, aes(x = year, y = SADsAvertedcum, colour = scenario)) +
    geom_line(linewidth = 0.6) +
    scale_y_continuous("Deaths avoided (cumulative)", labels = comma) +
    scale_x_continuous("Year", limits = c(2021, 2100),
                       breaks = seq(2030, 2100, 10)) +
    labs(title  = sprintf("%s — Smoke‑free Air‑Law Scenarios", state_name),
         colour = "Scenario") +
    theme_light(base_size = 11) +
    theme(legend.position = c(.05, .95),legend.justification = c("left", "top"), legend.box.just = "left",
          legend.background = element_rect(fill = alpha("white", 0.4), colour = NA),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 9),
          axis.text.x = element_text(angle = 60, hjust = 1))}

for (s in v_statefips) {
  sample_key <- grep(paste0("^", s, "_air_"), names(airlaws_out_list), value = TRUE)[1]
  abbr <- airlaws_out_list[[sample_key]][["df_mort.outputs"]]$abbr[1]
  state_name <- abbr_to_name(abbr)
  
  ggsave(filename = sprintf("figs/airlaws/%s_%s_DeathAvoided_airlaws.jpg", s, gsub(" ", "_", state_name)),
         plot = plot_DeathAvoided_airlaws(s), width = 6, height = 4, dpi = 600)}

# ------------------------------------------------------------------------------
# 3. Smoking Prevalence
# ------------------------------------------------------------------------------ 

plot_Prev_airlaws <- function(fipscode, scenarios = NULL, out_list = airlaws_out_list) {
  if (is.null(scenarios)) {
    scenarios <- apply(expand.grid(Iwp = 0:1, Ir = 0:1, Ib = 0:1), 1,
                       \(v) sprintf("wp%d_r%d_b%d", v[1], v[2], v[3]))
    scenarios <- scenarios[scenarios != "wp0_r0_b0"]}
  
  df_all <- lapply(scenarios, function(scen) {
    key <- sprintf("%s_air_%s", fipscode, scen)
    df <- out_list[[key]][["df_CSprevs.by.state"]]
    df$scenario <- scen
    df }) |> bind_rows()
  
  df_plot <- df_all |>
    filter(tolower(gender) == "both", year >= 2025) |>
    mutate(gender = "Persons") |>
    group_by(year, scenario) |>
    summarise(prev = mean(prev), .groups = "drop") |>
    group_by(scenario) |>
    filter(any(prev != 0) | scenario == "wp1_r1_b1") |>
    ungroup()
  
  scenario_levels <- c(
    "wp0_r0_b1", "wp0_r1_b0", "wp0_r1_b1",
    "wp1_r0_b0", "wp1_r0_b1", "wp1_r1_b0", "wp1_r1_b1")
  
  df_plot <- df_plot |>
    filter(scenario %in% scenario_levels) |>
    mutate(scenario = factor(scenario, levels = scenario_levels))
  
  state_name <- abbr_to_name(df_all$abbr[1])
  
  ggplot(df_plot, aes(x = year, y = prev, colour = scenario)) +
    geom_line(linewidth = 0.6) +
    scale_y_continuous("Smoking prevalence (%)",
                       labels = scales::percent_format(accuracy = 0.1)) +
    scale_x_continuous("Year", limits = c(2021, 2100),
                       breaks = seq(2030, 2100, 10)) +
    labs(title  = sprintf("%s — Smoking Prevalence", state_name),
         colour = "Scenario") +
    theme_light(base_size = 11) +
    theme(legend.position = c(.83, .95),  legend.justification = c("right", "top"),
          legend.background = element_rect(fill = alpha("white", .4), colour = NA),
          legend.title = element_text(size = 9), legend.text = element_text(size = 9),
          axis.text.x = element_text(angle = 60, hjust = 1))}

for (s in v_statefips) {
  sample_key <- grep(paste0("^", s, "_air_"), names(airlaws_out_list), value = TRUE)[1]
  abbr <- airlaws_out_list[[sample_key]][["df_mort.outputs"]]$abbr[1]
  state_name <- abbr_to_name(abbr)
  
  ggsave(filename = sprintf("figs/airlaws/%s_%s_Prev_airlaws.jpg", s, gsub(" ", "_", state_name)),
         plot = plot_Prev_airlaws(s), width = 6, height = 4, dpi = 600)}

# ------------------------------------------------------------------------------
# 4. Prevalence Reduction
# ------------------------------------------------------------------------------

normalise_gender <- function(x){
  tolower(x) |> recode("men" = "Males", "women" = "Females", "both" = "Persons")}

airlaws_prevDiff <- list()

scenarios  <- setdiff(
  apply(expand.grid(Iwp=0:1, Ir=0:1, Ib=0:1), 1,
        \(v) sprintf("wp%d_r%d_b%d", v[1], v[2], v[3])), "wp0_r0_b0")

for (s in v_statefips){
  base_tbl <- airlaws_out_list[[paste0(s, "_air_wp0_r0_b0")]][["df_CSprevs.by.state"]] %>%
    mutate(gender = normalise_gender(gender)) %>%
    select(year, gender, age, prev, abbr)
  
  for (sc in scenarios){
    pol_tbl <- airlaws_out_list[[paste0(s, "_air_", sc)]][["df_CSprevs.by.state"]] %>%
      mutate(gender = normalise_gender(gender)) %>%
      select(year, gender, age, prev)
    
    diff_tbl <- base_tbl %>%
      inner_join(pol_tbl, by = c("year","gender","age"),
                 suffix = c("_base","_pol")) %>%
      mutate(red_pp = (prev_base - prev_pol) * 100) %>%
      select(year, gender, age, red_pp, abbr)       
    
    airlaws_prevDiff[[paste0(s, "_diff_", sc)]] <- diff_tbl}}

plot_PrevRed_airlaws <- function(fipscode, scenarios = NULL, out_list = airlaws_prevDiff) {
  
  if (is.null(scenarios)) {
    scenarios <- apply(expand.grid(Iwp = 0:1, Ir = 0:1, Ib = 0:1), 1,
                       \(v) sprintf("wp%d_r%d_b%d", v[1], v[2], v[3]))
    scenarios <- scenarios[scenarios != "wp0_r0_b0"]}
  
  df_all <- lapply(scenarios, function(scen) {
    key <- sprintf("%s_diff_%s", fipscode, scen)
    df  <- out_list[[key]]
    df$scenario <- scen
    df }) |> bind_rows()
  
  df_plot <- df_all |>
    filter(gender == "Persons", year >= 2025) |>
    group_by(year, scenario) |>
    summarise(PrevRed = mean(red_pp), .groups = "drop") |>
    group_by(scenario) |>
    filter(any(PrevRed != 0) | scenario == "wp1_r1_b1") |>
    ungroup()
  
  scenario_levels <- c(
    "wp0_r0_b1", "wp0_r1_b0", "wp0_r1_b1",
    "wp1_r0_b0", "wp1_r0_b1", "wp1_r1_b0", "wp1_r1_b1")
  
  df_plot <- df_plot |>
    filter(scenario %in% scenario_levels) |>
    mutate(scenario = factor(scenario, levels = scenario_levels))
  
  abbr <- out_list[[sprintf("%s_diff_wp1_r1_b1", fipscode)]]$abbr[1]
  state_name <- abbr_to_name(abbr)
  
  ggplot(df_plot, aes(x = year, y = PrevRed, colour = scenario)) +
    geom_line(linewidth = 0.6) +
    scale_y_continuous("Prevalence reduction (pp)",
                       labels = \(x) sprintf("%.1f pp", x),
                       limits = c(0, NA)) +
    scale_x_continuous("Year", limits = c(2021, 2100),
                       breaks = seq(2030, 2100, 10)) +
    labs(title = sprintf("%s — Smoking Prevalence Reduction", state_name),
         colour = "Scenario",
         caption = "pp = percentage-points") +
    theme_light(base_size = 11) +
    theme(legend.position = "top", legend.direction  = "horizontal", legend.justification = "center", legend.box = "horizontal",     
          legend.margin = margin(b = 4),legend.background = element_rect(fill = alpha("white", 0.4), colour = NA),
          legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.spacing.x  = unit(6, "pt"),  
          legend.key.height = unit(9, "pt"), axis.text.x = element_text(angle = 60, hjust = 1))}

for (s in v_statefips) {
  abbr <- airlaws_prevDiff[[paste0(s, "_diff_wp1_r1_b1")]]$abbr[1]
  state_name <- abbr_to_name(abbr)
  
  ggsave(filename = sprintf("figs/airlaws/%s_%s_airlaw_PrevRed.jpg",s, gsub(" ", "_", state_name)),
         plot = plot_PrevRed_airlaws(s),
         width = 6, height = 4, dpi = 600)}






## ----------------------- Tobacco Control Expenditures Visualization ----------

# ------------------------------------------------------------------------------
# 1. Life Years Gained
# ------------------------------------------------------------------------------

plot_LYG_tcexp <- function(fipscode, out_list = tcexp_out_list) {
  scen_keys <- grep(paste0("^", fipscode, "_exp"), names(out_list), value = TRUE)
  
  df_all <- lapply(scen_keys, function(k) {
    df <- out_list[[k]][["df_mort.outputs"]]
    df$exp_scen <- sub("^\\d+_", "", k)  # keep "exp10", "exp20", …
    df }) |> bind_rows()
  
  df_plot <- df_all |>
    mutate(gender = recode(tolower(gender), "both" = "Persons")) |>
    filter(gender == "Persons", year >= 2025) |>
    group_by(exp_scen, year) |>
    summarise(LYGcum = sum(LYGcum), .groups = "drop")
  
  wide <- pivot_wider(df_plot, names_from = year, values_from = LYGcum)
  rep_scen <- wide |>
    group_by(across(-exp_scen)) |>
    summarise(chosen = first(exp_scen), .groups = "drop") |>
    pull(chosen)
  
  df_plot <- df_plot |>
    filter(exp_scen %in% rep_scen) |>
    mutate(exp_scen = factor(exp_scen, levels = rep_scen))
  
  state_name <- abbr_to_name(df_all$abbr[1])
  
  n_cols <- length(rep_scen)
  pal <- RColorBrewer::brewer.pal(max(3, min(8, n_cols)), "Dark2")
  
  ggplot(df_plot, aes(year, LYGcum, colour = exp_scen)) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual("Funding level\n(% CDC rec.)",
                        values = pal[seq_along(rep_scen)],
                        labels = sub("exp", "", rep_scen) |> paste0("%")) +
    scale_y_continuous("Life-years Gained", labels = scales::comma) +
    scale_x_continuous("Year", limits = c(2021, 2100),
                       breaks = seq(2030, 2100, 10)) +
    labs(title = sprintf("%s — Tobacco‑control spending scenarios", state_name)) +
    theme_light(base_size = 11) +
    theme(legend.position = c(.05, .95),legend.justification = c("left", "top"), legend.box.just = "left",
          legend.background = element_rect(fill = alpha("white", 0.4), colour = NA),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 9),
          legend.spacing.y = unit(6, "pt"), 
          legend.key.height = unit(9, "pt"),
          axis.text.x = element_text(angle = 60, hjust = 1))}

for (s in v_statefips) {
  scen_keys <- grep(paste0("^", s, "_exp"), names(tcexp_out_list), value = TRUE)
  state_name <- fips(as.numeric(s), to = "Name")
  ggsave(filename = sprintf("figs/tcexp/%s_%s_LYG_tcexp.jpg", s, gsub(" ", "_", state_name)),
         plot = plot_LYG_tcexp(s), width = 6, height = 4, dpi = 600)}

# ------------------------------------------------------------------------------
# 2. Death Avoided
# ------------------------------------------------------------------------------

plot_DeathAvoided_tcexp <- function(fipscode, out_list = tcexp_out_list) {
  scen_keys <- grep(paste0("^", fipscode, "_exp"), names(out_list), value = TRUE)
  
  df_all <- lapply(scen_keys, function(k) {
    df <- out_list[[k]][["df_mort.outputs"]]
    df$exp_scen <- sub("^\\d+_", "", k)  # keep "exp10", "exp20", …
    df }) |> bind_rows()
  
  df_plot <- df_all |>
    mutate(gender = recode(tolower(gender), "both" = "Persons")) |>
    filter(gender == "Persons", year >= 2025) |>
    group_by(exp_scen, year) |>
    summarise(SADsAvertedcum = sum(SADsAvertedcum), .groups = "drop")
  
  wide <- pivot_wider(df_plot, names_from = year, values_from = SADsAvertedcum)
  rep_scen <- wide |>
    group_by(across(-exp_scen)) |>
    summarise(chosen = first(exp_scen), .groups = "drop") |>
    pull(chosen)
  
  df_plot <- df_plot |>
    filter(exp_scen %in% rep_scen) |>
    mutate(exp_scen = factor(exp_scen, levels = rep_scen))
  
  state_name <- abbr_to_name(df_all$abbr[1])
  
  n_cols <- length(rep_scen)
  pal <- RColorBrewer::brewer.pal(max(3, min(8, n_cols)), "Dark2")
  
  ggplot(df_plot, aes(year,SADsAvertedcum, colour = exp_scen)) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual("Funding level\n(% CDC rec.)",
                        values = pal[seq_along(rep_scen)],
                        labels = sub("exp", "", rep_scen) |> paste0("%")) +
    scale_y_continuous("Deaths avoided (cumulative)", labels = comma) +
    scale_x_continuous("Year", limits = c(2021, 2100),
                       breaks = seq(2030, 2100, 10)) +
    labs(title = sprintf("%s — Tobacco‑control spending scenarios", state_name)) +
    theme_light(base_size = 11) +
    theme(legend.position = c(.05, .95),legend.justification = c("left", "top"), legend.box.just = "left",
          legend.background = element_rect(fill = alpha("white", 0.4), colour = NA),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 9),
          legend.spacing.y = unit(6, "pt"), 
          legend.key.height = unit(9, "pt"),
          axis.text.x = element_text(angle = 60, hjust = 1))}

for (s in v_statefips) {
  scen_keys <- grep(paste0("^", s, "_exp"), names(tcexp_out_list), value = TRUE)
  state_name <- fips(as.numeric(s), to = "Name")
  ggsave(filename = sprintf("figs/tcexp/%s_%s_DeathsAvoided_tcexp.jpg", s, gsub(" ", "_", state_name)),
         plot = plot_DeathAvoided_tcexp(s), width = 6, height = 4, dpi = 600)}


# ------------------------------------------------------------------------------
# 3. Smoking Prevalence
# ------------------------------------------------------------------------------ 

plot_Prev_tcexp <- function(fipscode, out_list = tcexp_out_list) {
  scen_keys <- grep(paste0("^", fipscode, "_exp"), names(out_list), value = TRUE)
  
  df_all <- lapply(scen_keys, function(k) {
    df <- out_list[[k]][["df_CSprevs.by.state"]]
    df$exp_scen <- sub("^\\d+_", "", k)  # keep "exp10", "exp20", …
    df }) |> bind_rows()
  
  df_plot <- df_all |>
    mutate(gender = recode(tolower(gender), "both" = "Persons")) |>
    filter(gender == "Persons", year >= 2025) |>
    group_by(exp_scen, year) |>
    summarise(prev_pct = mean(prev) * 100, .groups = "drop")
  
  wide <- pivot_wider(df_plot, names_from = year, values_from = prev_pct)
  rep_scen <- wide |>
    group_by(across(-exp_scen)) |>
    summarise(chosen = first(exp_scen), .groups = "drop") |>
    pull(chosen)
  
  df_plot <- df_plot |>
    filter(exp_scen %in% rep_scen) |>
    mutate(exp_scen = factor(exp_scen, levels = rep_scen))
  
  state_name <- abbr_to_name(df_all$abbr[1])
  
  n_cols <- length(rep_scen)
  pal <- RColorBrewer::brewer.pal(max(3, min(8, n_cols)), "Dark2")
  
  ggplot(df_plot, aes(year,prev_pct, colour = exp_scen)) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual("Funding level\n(% CDC rec.)",
                        values = pal[seq_along(rep_scen)],
                        labels = sub("exp", "", rep_scen) |> paste0("%")) +
    scale_y_continuous("Smoking prevalence (%)",
                       labels = function(x) sprintf("%.1f%%", x)) +
    scale_x_continuous("Year", limits = c(2021, 2100),
                       breaks = seq(2030, 2100, 10)) +
    labs(title = sprintf("%s — Tobacco‑control spending scenarios", state_name)) +
    theme_light(base_size = 11) +
    theme(legend.position = c(.85, .95),legend.justification = c("right", "top"), legend.box.just = "left",
          legend.background = element_rect(fill = alpha("white", 0.4), colour = NA),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 6),
          legend.spacing.y = unit(6, "pt"), 
          legend.key.height = unit(6, "pt"),
          axis.text.x = element_text(angle = 60, hjust = 1))}

for (s in v_statefips) {
  scen_keys <- grep(paste0("^", s, "_exp"), names(tcexp_out_list), value = TRUE)
  state_name <- fips(as.numeric(s), to = "Name")
  ggsave(filename = sprintf("figs/tcexp/%s_%s_Prev_tcexp.jpg", s, gsub(" ", "_", state_name)),
         plot = plot_Prev_tcexp(s), width = 6, height = 4, dpi = 600)}


# ------------------------------------------------------------------------------
# 4. Prevalence Reduction
# ------------------------------------------------------------------------------

normalise_gender <- function(x){
  tolower(x) |> recode("male" = "Males", "female" = "Females","both" = "Persons")}

tcexp_prevDiff <- list()

for (s in v_statefips){
  base_tbl <- baseline_results[[s]][["df_CSprevs.by.state"]] %>%
    mutate(gender = normalise_gender(gender)) %>%
    select(year, gender, age, prev_base = prev)
  
  scen_keys <- grep(paste0("^", s, "_exp"), names(tcexp_out_list), value = TRUE)
  
  for (key in scen_keys){
    
    pol_tbl <- tcexp_out_list[[key]][["df_CSprevs.by.state"]] %>%
      mutate(gender = normalise_gender(gender)) %>%
      select(year, gender, age, prev_tax = prev)
    
    diff_tbl <- base_tbl %>%
      inner_join(pol_tbl, by = c("year","gender","age")) %>%
      mutate(prev_diff_pp = (prev_base - prev_tax) * 100) %>%
      select(year, gender, age, prev_diff_pp)
    
    tcexp_prevDiff[[paste0(key, "_diff")]] <- diff_tbl   # e.g. 01_exp30_diff
  }
}

plot_PrevRed_tcexp <- function(fipscode, out_list_diff = tcexp_prevDiff) {
  scen_keys <- grep(paste0("^", fipscode, "_exp.*_diff$"), names(out_list_diff), value = TRUE)
  if (length(scen_keys) == 0) return(NULL)
  
  df_all <- lapply(scen_keys, function(k) {
    df <- out_list_diff[[k]]
    df$exp_scen <- sub("^\\d+_", "", k)  # keep "exp10", etc.
    df
  }) |> bind_rows()
  
  df_plot <- df_all |>
    filter(tolower(gender) == "persons", year >= 2025) |>
    group_by(exp_scen, year) |>
    summarise(prev_red_pp = mean(prev_diff_pp), .groups = "drop")
  
  wide <- pivot_wider(df_plot, names_from = year, values_from = prev_red_pp)
  rep_scen <- wide |>
    group_by(across(-exp_scen)) |>
    summarise(chosen = first(exp_scen), .groups = "drop") |>
    pull(chosen)
  
  df_plot <- df_plot |>
    filter(exp_scen %in% rep_scen) |>
    mutate(exp_scen = factor(exp_scen, levels = rep_scen))
  
  # Retrieve state name from any baseline source
  abbr <- baseline_results[[fipscode]][["df_CSprevs.by.state"]]$abbr[1]
  state_name <- abbr_to_name(abbr)
  
  n_cols <- length(rep_scen)
  pal <- RColorBrewer::brewer.pal(max(3, min(8, n_cols)), "Dark2")
  
  ggplot(df_plot, aes(year, prev_red_pp, colour = exp_scen)) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual("Funding level\n(% CDC rec.)",
                        values = pal[seq_along(rep_scen)],
                        labels = gsub("(exp|_diff)", "", rep_scen) |> paste0("%")) +
    scale_y_continuous("Prevalence reduction (pp)",
                       labels = \(x) sprintf("%.1f pp", x)) +
    scale_x_continuous("Year", limits = c(2021, 2100),
                       breaks = seq(2030, 2100, 10)) +
    labs(title = sprintf("%s — Tobacco‑control spending scenarios", state_name),
         caption = "pp = percentage-points") +
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

for (s in v_statefips) {
  state_name <- fips(as.numeric(s), to = "Name")
  ggsave(filename = sprintf("figs/tcexp/%s_%s_PrevRed_tcexp.jpg", s, gsub(" ", "_", state_name)),
         plot = plot_PrevRed_tcexp(s), width = 6, height = 4, dpi = 600)}




## ------------ Health Graphic Warnings (policy year = 2025) Visualization -----

# ------------------------------------------------------------------------------
# 1. Life Years Gained
# ------------------------------------------------------------------------------

plot_LYG_warnings <- function(fipscode, warnings_out_list_2025) {
  scen_keys <- grep(paste0("^", fipscode, "_warn"), names(warnings_out_list_2025), value = TRUE)
  
  df_all <- lapply(scen_keys, function(k) {
    df <- warnings_out_list_2025[[k]][["df_mort.outputs"]]
    df$scen <- k
    df }) |> bind_rows()
  
  df_plot <- df_all |> 
    mutate(gender = tolower(gender)) |> 
    filter(gender == "both", year >= 2025) |> 
    group_by(scen, year) |> 
    summarise(LYGcum = sum(LYGcum), .groups = "drop")
  
  top_scen <- df_plot |> 
    filter(year == 2100) |> 
    arrange(desc(LYGcum)) |> 
    slice_head(n = 3) |> 
    pull(scen)
  
  df_plot <- df_plot |> 
    filter(scen %in% top_scen) |> 
    mutate(scen = factor(scen, levels = top_scen))
  
  abbr <- df_all$abbr[1]
  state_name <- abbr_to_name(abbr)
  
  n_cols <- length(top_scen)
  pal <- RColorBrewer::brewer.pal(max(3, min(8, n_cols)), "Set1")
  
  ggplot(df_plot, aes(year, LYGcum, colour = scen)) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual("Scenario",
                        values = pal[seq_along(top_scen)],
                        labels = gsub("^\\d+_warn_", "", top_scen)) +
    scale_y_continuous("Life-years Gained", labels = scales::comma) +
    scale_x_continuous("Year", limits = c(2021, 2100), breaks = seq(2030, 2100, 10)) +
    labs(title = sprintf("%s — Graphic Health Warning (Top 3 Scenarios)", state_name)) +
    theme_light(base_size = 11) +
    theme(legend.position = c(0.05, 0.90),legend.justification = c("left", "top"),legend.box.just = "left",
          legend.background = element_rect(fill = alpha("white", 0.5), colour = NA),
          legend.title = element_text(size = 9), 
          legend.text = element_text(size = 8),
          axis.text.x = element_text(angle = 60, hjust = 1))}

for (s in v_statefips) {
  state_name <- fips(as.numeric(s), to = "Name")
  ggsave(filename = sprintf("figs/warnings/%s_%s_LYG_warnings.jpg", s, gsub(" ", "_", state_name)),
         plot = plot_LYG_warnings(s, warnings_out_list_2025), width = 6, height = 4, dpi = 600)}


# ------------------------------------------------------------------------------
# 2. Death Avoided
# ------------------------------------------------------------------------------

plot_DeathAvoided_warnings <- function(fipscode, warnings_out_list_2025) {
  scen_keys <- grep(paste0("^", fipscode, "_warn"), names(warnings_out_list_2025), value = TRUE)
  
  df_all <- lapply(scen_keys, function(k) {
    df <- warnings_out_list_2025[[k]][["df_mort.outputs"]]
    df$scen <- k
    df }) |> bind_rows()
  
  df_plot <- df_all |> 
    mutate(gender = tolower(gender)) |> 
    filter(gender == "both", year >= 2025) |> 
    group_by(scen, year) |> 
    summarise(SADsAvertedcum = sum(SADsAvertedcum), .groups = "drop")
  
  top_scen <- df_plot |> 
    filter(year == 2100) |> 
    arrange(desc(SADsAvertedcum)) |> 
    slice_head(n = 3) |> 
    pull(scen)
  
  df_plot <- df_plot |> 
    filter(scen %in% top_scen) |> 
    mutate(scen = factor(scen, levels = top_scen))
  
  abbr <- df_all$abbr[1]
  state_name <- abbr_to_name(abbr)
  
  n_cols <- length(top_scen)
  pal <- RColorBrewer::brewer.pal(max(3, min(8, n_cols)), "Set1")
  
  ggplot(df_plot, aes(year, SADsAvertedcum, colour = scen)) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual("Scenario",
                        values = pal[seq_along(top_scen)],
                        labels = gsub("^\\d+_warn_", "", top_scen)) +
    scale_y_continuous("Deaths avoided (cumulative)", labels = scales::comma) +
    scale_x_continuous("Year", limits = c(2021, 2100), breaks = seq(2030, 2100, 10)) +
    labs(title = sprintf("%s — Graphic Health Warning (Top 3 Scenarios)", state_name)) +
    theme_light(base_size = 11) +
    theme(legend.position = c(0.05, 0.90),legend.justification = c("left", "top"),legend.box.just = "left",
          legend.background = element_rect(fill = alpha("white", 0.5), colour = NA),
          legend.title = element_text(size = 9), 
          legend.text = element_text(size = 8),
          axis.text.x = element_text(angle = 60, hjust = 1))}

for (s in v_statefips) {
  state_name <- fips(as.numeric(s), to = "Name")
  ggsave(filename = sprintf("figs/warnings/%s_%s_DeathsAvoided_warnings.jpg", s, gsub(" ", "_", state_name)),
         plot = plot_DeathAvoided_warnings(s, warnings_out_list_2025), width = 6, height = 4, dpi = 600)}


# ------------------------------------------------------------------------------
# 3. Smoking Prevalence
# ------------------------------------------------------------------------------ 

plot_Prev_warnings <- function(fipscode, warnings_out_list_2025) {
  scen_keys <- grep(paste0("^", fipscode, "_warn"), names(warnings_out_list_2025), value = TRUE)
  
  df_all <- lapply(scen_keys, function(k) {
    df <- warnings_out_list_2025[[k]][["df_CSprevs.by.state"]]
    df$scen <- k
    df }) |> bind_rows()
  
  df_plot <- df_all |> 
    mutate(gender = tolower(gender)) |> 
    filter(gender == "both", year >= 2025) |> 
    group_by(scen, year) |> 
    summarise(prev_pct = mean(prev) * 100, .groups = "drop")
  
  top_scen <- df_plot |> 
    filter(year == 2100) |> 
    arrange(prev_pct) |> 
    slice_head(n = 3) |> 
    pull(scen)
  
  df_plot <- df_plot |> 
    filter(scen %in% top_scen) |> 
    mutate(scen = factor(scen, levels = top_scen))
  
  abbr <- df_all$abbr[1]
  state_name <- abbr_to_name(abbr)
  
  n_cols <- length(top_scen)
  pal <- RColorBrewer::brewer.pal(max(3, min(8, n_cols)), "Set1")
  
  ggplot(df_plot, aes(year, prev_pct, colour = scen)) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual("Scenario",
                        values = pal[seq_along(top_scen)],
                        labels = gsub("^\\d+_warn_", "", top_scen)) +
    scale_y_continuous("Smoking prevalence (%)", labels = function(x) sprintf("%.1f%%", x)) +
    scale_x_continuous("Year", limits = c(2021, 2100),
                       breaks = seq(2030, 2100, 10)) +
    labs(title = sprintf("%s — Graphic Health Warning (Top 3 Scenarios)", state_name),
         caption = "Persons = Both sexes") +
    theme_light(base_size = 11) +
    theme(legend.position = c(.65, .95),legend.justification = c("left", "top"), legend.box.just = "left",
          legend.background = element_rect(fill = alpha("white", 0.4), colour = NA),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 9),
          axis.text.x = element_text(angle = 60, hjust = 1))}

for (s in v_statefips) {
  state_name <- fips(as.numeric(s), to = "Name")
  ggsave(filename = sprintf("figs/warnings/%s_%s_Prev_warnings.jpg", s, gsub(" ", "_", state_name)),
         plot = plot_Prev_warnings(s, warnings_out_list_2025), width = 6, height = 4, dpi = 600)}


# ------------------------------------------------------------------------------
# 4. Prevalence Reduction
# ------------------------------------------------------------------------------

normalise_gender <- function(x){
  tolower(x) |> recode("men" = "Males", "women" = "Females","both" = "Persons")}

warnings_prevDiff <- list()

for (s in v_statefips){
  base_tbl <- baseline_results[[s]][["df_CSprevs.by.state"]] %>%
    mutate(gender = normalise_gender(gender)) %>%
    select(year, gender, age, prev_base = prev)
  
  scen_keys <- grep(paste0("^", s, "_warn_"), names(warnings_out_list_2025), value = TRUE)
  
  
  for (key in scen_keys){
    pol_tbl <- warnings_out_list_2025[[key]][["df_CSprevs.by.state"]] %>%
      mutate(gender = normalise_gender(gender)) %>%
      select(year, gender, age, prev_tax = prev)
    
    diff_tbl <- base_tbl %>%
      inner_join(pol_tbl, by = c("year","gender","age")) %>%
      mutate(prev_diff_pp = (prev_base - prev_tax) * 100) %>%
      select(year, gender, age, prev_diff_pp)
    warnings_prevDiff[[paste0(key, "_diff")]] <- diff_tbl 
  }
}

plot_PrevRed_warnings <- function(fipscode, out_list_diff = warnings_prevDiff) {
  scen_keys <- grep(paste0("^", fipscode, "_warn_.*_diff$"), names(out_list_diff), value = TRUE)
  if (length(scen_keys) == 0) return(NULL)
  
  df_all <- lapply(scen_keys, function(k) {
    df <- out_list_diff[[k]]
    df$scen <- k |> sub("^\\d+_warn_", "", x = _) |> sub("_diff$", "", x = _)
    df }) |> bind_rows()
  
  df_plot <- df_all |>
    filter(gender == "Persons", year >= 2022) |>
    group_by(scen, year) |>
    summarise(prev_pct = mean(prev_diff_pp), .groups = "drop")
  
  top_scen <- df_plot |>
    filter(year == 2100) |>
    arrange(desc(prev_pct)) |>  # higher prev_diff_pp = better
    slice_head(n = 3) |>
    pull(scen)
  
  df_plot <- df_plot |>
    filter(scen %in% top_scen) |>
    mutate(scen = factor(scen, levels = top_scen))
  
  abbr <- baseline_results[[fipscode]][["df_CSprevs.by.state"]]$abbr[1]
  state_name <- abbr_to_name(abbr)
  
  pal <- RColorBrewer::brewer.pal(max(3, length(top_scen)), "Dark2")
  
  ggplot(df_plot, aes(year, prev_pct, colour = scen)) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual("Scenario", values = pal[seq_along(top_scen)]) +
    scale_y_continuous("Prevalence reduction (pp)", labels = \(x) sprintf("%.1f pp", x)) +
    scale_x_continuous("Year", limits = c(2021, 2100), breaks = seq(2030, 2100, 10)) +
    labs(title = sprintf("%s — Graphic Health Warning (Top 3 Scenarios)", state_name),
         caption = "pp = percentage-points") +
    theme_light(base_size = 11) +
    theme(legend.position = c(.65, .65),legend.justification = c("left", "top"), legend.box.just = "left",
          legend.background = element_rect(fill = alpha("white", 0.4), colour = NA),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 9),
          axis.text.x = element_text(angle = 60, hjust = 1))}

for (s in v_statefips) {
  state_name <- fips(as.numeric(s), to = "Name")
  ggsave(filename = sprintf("figs/warnings/%s_%s_PrevRed_warnings.jpg", s, gsub(" ", "_", state_name)),
         plot = plot_PrevRed_warnings(s, warnings_prevDiff), width = 6, height = 4, dpi = 600)}










