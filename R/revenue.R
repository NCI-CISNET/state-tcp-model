

## =============================================================================
## Tax Revenue Analysis Based on CoRRE
## =============================================================================

library(dplyr)
library(readr)
library(tidyr) 
library(tidycensus)
library(purrr)
library(readxl)
library(tidycensus)

# =============================================================================
## Parameters
date_variable <- format(Sys.Date(), "%m.%d.%y")
base_year <- 2025
policy_year <- 2026
price_elasticity <- -0.2 

## Load Model Results ----
# baseline_results <- readRDS("./analysis_output/20250510/baseline_results.rds") # TCP Model Baseline results
# tax_out_list <- readRDS("./analysis_output/20250510/tax_out_list.rds") # TCP Model Tax results (new elasticity, policyyear=2026)

## Load External Data Sources
# Pack Sales Per Capita; Year = 2019 (Cigarette consumption)
# packs_pc <- read_xlsx("./tax_prevalence_economics/cigarette_consumption_pack_sales_per_capita.xlsx") 
# packs_pc  <- packs_pc %>% filter(Year == 2019)

# Baseline cigarette price and tax
# baseline_price_per_pack <- read.csv("./tax_prevalence_economics/cigarette_prices_by_state_2025.csv") # baseline price per pack: Year 2025
# tax_per_pack <- read.csv("./tax_prevalence_economics/cigarette_prices_by_state_2025.csv") # state tax per pack Year: Year 2025



baseline_results <- readRDS("./R/baseline_results.rds") # TCP Model Baseline results
tax_out_list <- readRDS("./R/tax_out_list.rds")# TCP Model Tax results (new elasticity, policyyear=2026)
packs_pc <- read_xlsx("./R/cigarette_consumption_pack_sales_per_capita.xlsx") 
packs_pc  <- packs_pc %>% filter(Year == 2019)
baseline_price_per_pack <- read.csv("./R/cigarette_prices_by_state_2025.csv") # baseline price per pack: Year 2025
state_tax_per_pack <- read_xlsx("./R/state_tax_per_pack.xlsx") # state tax per pack: latest available year by state


# =============================================================================
# State population data PROCESSING: From census population data; Year: 2025 - 2100

state_pop <- list() 
for (f in v_statefips) {
  load(file = paste0("./data/state_inputs/pop_", f, ".RData"))  # load df_F.census_data + df_M.census_data
  state_pop[[f]] <- list(female = df_F.census_data, male = df_M.census_data)}

## Calculate Total State Population by Year ----
state_total_pop <- lapply(names(state_pop), function(fips) { 
  df_female <- state_pop[[fips]]$female
  df_male <- state_pop[[fips]]$male
  
  female_total_by_year <- colSums(df_female, na.rm = TRUE)
  male_total_by_year <- colSums(df_male, na.rm = TRUE)
  
  # Convert from wide to long format
  df_female_long <- df_female %>% pivot_longer(cols = everything(), names_to = "year", values_to = "female") %>% mutate(year = as.integer(year))
  df_male_long <- df_male %>% pivot_longer(cols = everything(), names_to = "year", values_to = "male") %>% mutate(year = as.integer(year))
  
  df_total <- tibble(
    year = as.integer(names(female_total_by_year)),
    female = as.numeric(female_total_by_year),
    male = as.numeric(male_total_by_year),
    fips = fips
  ) %>%
    mutate(total_population = female + male) %>%
    select(year, fips, total_population) %>%
    filter(year >= base_year)              # <-- keep 2025+
}) %>% bind_rows()

## Add State Identifiers ----
state_lookup <- fips_codes %>%
  select(state_name, state_code, state) %>%
  distinct() %>%
  rename(abbr = state, fips = state_code)

state_total_pop <- state_total_pop %>% left_join(state_lookup,by = c("fips"))

# =============================================================================
# Data imputation

## Time Range (Year2025 - 2100)
years_full <- seq(base_year, max(state_total_pop$year))

## Process Cigarette Consumption Data (Year2019) ---- 
packs_clean <- packs_pc %>%
  transmute(
    year = as.integer(Year),
    abbr = LocationAbbr,
    state_name = LocationDesc,                      # from packs
    pack_sales_per_capita = as.numeric(Data_Value)
  ) %>%
  left_join(state_lookup, by = "abbr", suffix = c("", "_lkp")) %>%  # avoid silent suffixing
  mutate(
    state_name = coalesce(state_name, state_name_lkp)
  ) %>%
  select(year, abbr, fips, state_name, pack_sales_per_capita) %>%
  group_by(fips, abbr) %>%
  complete(year = years_full) %>%
  arrange(year) %>%
  fill(state_name, .direction = "downup") %>%       # label fill both ways
  fill(pack_sales_per_capita, .direction = "down") %>%  # LOCF to 2100
  ungroup()

## Process Tax and Price Data ----
lk <- setNames(state_lookup$abbr, state_lookup$state_name)
baseline_price_per_pack <- baseline_price_per_pack %>%
  mutate(State = dplyr::coalesce(unname(lk[trimws(State)]), toupper(trimws(State))))

price_clean <- baseline_price_per_pack %>%
  rename(abbr = State) %>%
  mutate(
    year = 2025L,
    baseline_price_per_pack = parse_number(as.character(price_per_pack_incl_taxes))
  ) %>%
  left_join(state_lookup, by = "abbr") %>%       # adds fips
  select(year, fips, abbr, baseline_price_per_pack) %>%
  group_by(fips, abbr) %>%
  complete(year = years_full) %>%
  arrange(year) %>%
  fill(baseline_price_per_pack, .direction = "down") %>%
  ungroup()

tax_clean <- state_tax_per_pack %>%
  transmute(
    abbr = trimws(LocationAbbr),
    year = as.integer(Year),
    tax_per_pack = as.numeric(Data_Value)
  ) %>%
  filter(!is.na(abbr), !is.na(year), !is.na(tax_per_pack)) %>%
  group_by(abbr) %>%
  slice_max(order_by = year, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(year = 2025L) %>%
  left_join(state_lookup, by = "abbr") %>%
  select(year, fips, abbr, tax_per_pack) %>%
  group_by(fips, abbr) %>%
  complete(year = years_full) %>%
  arrange(year) %>%
  fill(tax_per_pack, .direction = "down") %>%
  ungroup()

# =============================================================================
# ---- Baseline (consumption, revenue, intensity) ----
# =============================================================================

## Calculate Baseline Consumption and Revenue ----
df_baseline <- state_total_pop %>%
  left_join(packs_clean, by = c("fips","abbr","year","state_name")) %>%
  left_join(tax_clean, by = c("fips","abbr","year")) %>%
  left_join(price_clean, by = c("fips","abbr","year")) %>%
  mutate(
    # Baseline total cigarette sales
    total_cigarette_consumption = total_population * pack_sales_per_capita, 
    # Baseline cigarette tax revenue
    baseline_state_tax_revenue  = total_cigarette_consumption * tax_per_pack) 

## Calculate Baseline Smoking Prevalence and Intensity ----
baseline_smoker <- purrr::map_dfr(names(baseline_results), function(fips_chr){
  
  # Census population
  census_pop_male <- state_pop[[fips_chr]]$male
  census_pop_female <- state_pop[[fips_chr]]$female
  
  # Smoking Prevalence from TCP Model
  l <- baseline_results[[fips_chr]]$l_prev_out
  cy <- readr::parse_number(colnames(l$m_M_popAP)) %>% as.integer()
  idx <- match(years_full, cy)
  
  # TCP Model population data to get smoking rate
  sm_model  <- (l$m_M_smokers + l$m_F_smokers)[, idx, drop = FALSE]
  pop_model <- (l$m_M_popAP + l$m_F_popAP  )[, idx, drop = FALSE]
  smoking_rate <- colSums(sm_model, na.rm = TRUE) / colSums(pop_model, na.rm = TRUE)
  
  census_years <- as.character(years_full)
  available_years <- intersect(census_years, colnames(census_pop_male))
  
  if(length(available_years) > 0) {
    census_total_pop <- colSums(census_pop_male[, available_years, drop = FALSE], na.rm = TRUE) + 
      colSums(census_pop_female[, available_years, drop = FALSE], na.rm = TRUE)
    
    census_total_pop <- census_total_pop[match(census_years, names(census_total_pop))]
  } else {
   
    census_total_pop <- colSums(pop_model, na.rm = TRUE)
  }
  
  tibble(fips = fips_chr, year = years_full,
         baseline_num_smokers = census_total_pop * smoking_rate,
         baseline_prev = smoking_rate,
         census_total_pop = census_total_pop)
})

## Merge Baseline Data and Calculate Smoking Intensity ----
df_baseline <- df_baseline %>%
  left_join(baseline_smoker, by = c("fips","year")) %>%
  mutate(
    # Baseline smoking intensity = Baseline total cigarette sales / Baseline number of adult smokers
    baseline_avg_smoking_intensity = ifelse(baseline_num_smokers > 0,
                                     total_cigarette_consumption / baseline_num_smokers,
                                     NA_real_))

# =============================================================================
#  POLICY SCENARIO CALCULATION - tax=3
# =============================================================================
tax <- 3

## Calculate Price Increase Effects ----
# Required percentage increase in cigarette price
df_policy <- df_baseline %>%
  mutate(pct_increase_price_per_pack = if_else(year >= policy_year & baseline_price_per_pack > 0,
                                               tax / baseline_price_per_pack, 0),
         pct_reduction_intensity = price_elasticity * pct_increase_price_per_pack)

## Extract Policy Smoking Prevalence ----
# Expected reduced number of smokers
policy_prev_tbl <- purrr::map_dfr(
  names(tax_out_list)[grepl("_3dollar$", names(tax_out_list))],
  function(tag){
    l  <- tax_out_list[[tag]]$l_prev_out
    cy <- readr::parse_number(colnames(l$m_M_popAP)) %>% as.integer()
    idx <- match(years_full, cy)
    
    ## From TCP Model to get $3tax increase smoking prevalence
    sm  <- (l$m_M_smokers + l$m_F_smokers)[, idx, drop = FALSE]
    pop <- (l$m_M_popAP   + l$m_F_popAP  )[, idx, drop = FALSE]
    tibble(fips = sub("_.*$","",tag), year = years_full,
           policy_prev = colSums(sm, na.rm = TRUE) / colSums(pop, na.rm = TRUE))
  })


## Calculate Final Policy Effects and Revenue ----
df_final_effects <- df_baseline %>%
  left_join(policy_prev_tbl, by = c("fips","year")) %>%
  #left_join(baseline_smoker %>% select(fips, year, census_total_pop), by = c("fips","year")) %>%  
  mutate(
    # Price increase percentage
    pct_increase_price_per_pack = if_else(year >= policy_year & baseline_price_per_pack > 0,
                                          tax / baseline_price_per_pack, 0),
    
    # Smoking intensity reduction (CoRRE elasticity: -0.20)
    pct_reduction_intensity = price_elasticity * pct_increase_price_per_pack,
    
    # Policy smoking intensity = Baseline intensity × (1 + % reduction in intensity)
    policy_smoking_intensity = baseline_avg_smoking_intensity * (1 + pct_reduction_intensity),
    
    # Policy number of smokers = Policy prevalence × Census Population
    policy_num_smokers = if_else(year >= policy_year,
                                 policy_prev * census_total_pop,  
                                 baseline_num_smokers),
    
    # Expected reduced annual total cigarette sales = Policy smokers × Policy intensity
    policy_total_cig_consumption = if_else(year >= policy_year,
                                           policy_num_smokers * policy_smoking_intensity,
                                           total_cigarette_consumption),
    
    # Expected cigarette tax revenue
    policy_state_tax_revenue = policy_total_cig_consumption *
      if_else(year >= policy_year, tax_per_pack + tax, tax_per_pack),
    
    # Revenue gain
    revenue_gain = if_else(year >= policy_year,
                           policy_state_tax_revenue - baseline_state_tax_revenue, 0)
  )

write.csv(df_final_effects, paste0("./R/$3_revenue_output_", date_variable, ".csv"), row.names = FALSE)


# =============================================================================
#  POLICY SCENARIO CALCULATION - tax=2
# =============================================================================
tax <- 2

## Calculate Price Increase Effects ----
# Required percentage increase in cigarette price
df_policy <- df_baseline %>%
  mutate(pct_increase_price_per_pack = if_else(year >= policy_year & baseline_price_per_pack > 0,
                                               tax / baseline_price_per_pack, 0),
         pct_reduction_intensity = price_elasticity * pct_increase_price_per_pack)

## Extract Policy Smoking Prevalence ----
# Expected reduced number of smokers
policy_prev_tbl <- purrr::map_dfr(
  names(tax_out_list)[grepl("_2dollar$", names(tax_out_list))],
  function(tag){
    l  <- tax_out_list[[tag]]$l_prev_out
    cy <- readr::parse_number(colnames(l$m_M_popAP)) %>% as.integer()
    idx <- match(years_full, cy)
    
    ## From TCP Model to get $3tax increase smoking prevalence
    sm  <- (l$m_M_smokers + l$m_F_smokers)[, idx, drop = FALSE]
    pop <- (l$m_M_popAP   + l$m_F_popAP  )[, idx, drop = FALSE]
    tibble(fips = sub("_.*$","",tag), year = years_full,
           policy_prev = colSums(sm, na.rm = TRUE) / colSums(pop, na.rm = TRUE))
  })


## Calculate Final Policy Effects and Revenue ----
df_final_effects <- df_baseline %>%
  left_join(policy_prev_tbl, by = c("fips","year")) %>%
  #left_join(baseline_smoker %>% select(fips, year, census_total_pop), by = c("fips","year")) %>%  
  mutate(
    # Price increase percentage
    pct_increase_price_per_pack = if_else(year >= policy_year & baseline_price_per_pack > 0,
                                          tax / baseline_price_per_pack, 0),
    
    # Smoking intensity reduction (CoRRE elasticity: -0.20)
    pct_reduction_intensity = price_elasticity * pct_increase_price_per_pack,
    
    # Policy smoking intensity = Baseline intensity × (1 + % reduction in intensity)
    policy_smoking_intensity = baseline_avg_smoking_intensity * (1 + pct_reduction_intensity),
    
    # Policy number of smokers = Policy prevalence × Census Population
    policy_num_smokers = if_else(year >= policy_year,
                                 policy_prev * census_total_pop,  
                                 baseline_num_smokers),
    
    # Expected reduced annual total cigarette sales = Policy smokers × Policy intensity
    policy_total_cig_consumption = if_else(year >= policy_year,
                                           policy_num_smokers * policy_smoking_intensity,
                                           total_cigarette_consumption),
    
    # Expected cigarette tax revenue
    policy_state_tax_revenue = policy_total_cig_consumption *
      if_else(year >= policy_year, tax_per_pack + tax, tax_per_pack),
    
    # Revenue gain
    revenue_gain = if_else(year >= policy_year,
                           policy_state_tax_revenue - baseline_state_tax_revenue, 0)
  )

write.csv(df_final_effects, paste0("./R/$2_revenue_output_", date_variable, ".csv"), row.names = FALSE)























