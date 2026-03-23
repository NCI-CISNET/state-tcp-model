### CISNET Population Model

rm(list = ls()) 
setwd("/Users/wangmengyao/Desktop/GitHub/state-tcp-model/")

# Package names
packages <- c("ggplot2", "readxl", "readr", "dplyr", "tidyr", "reshape2","grid", 
              "ggpubr", "gridBase", "gridExtra", "cdlTools", "stringr", "ggrepel")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
lapply(packages, library, character.only = TRUE)

#---------------------- Load data sets --------------------------------------------------------

for (f in v_statefips) {
  load(paste0("data/state_inputs/le_",f,".RData"))
  load(paste0("data/state_inputs/pop_",f,".RData"))
  load(paste0("data/state_inputs/mort_rates/mort_",f,".RData"))
  load(paste0("data/state_inputs/mort_rates/p.mort_", f, ".RData"))
  load(paste0("data/state_inputs/smk_",f,".RData"))
}

#-------------------- Set model inputs -----------------------------------------------------------
##### Scenario indicators 
##### policy_decay -> 0/1 default is standard policy effects (no exponential decay) =0 
##### make_tcp_out -> 0/1 default is not to make TCP tool files =0
policy_decay <- 0
make_tcp_out <- 0

#other model inputs

date_variable <- format(Sys.Date(), "%m.%d.%y")
startbc <- 1908   # starting birth cohort 
endbc <- 2100     # ending birth cohort
endyear <- 2200   # final calendar year, extra 100 years needed for AC-AP conversion

cohyears <- endbc-startbc+1     # number of cohort years
totalyears<-length(startbc:endyear)
v_calyears <- (startbc-startbc+1):(endbc-startbc+1)  # index of calendar years 1908-2100
v_stdbirths<- rep(1000000, times = 193)    #fixed birthrate/population size for prevalence calculations

# vector of states to run simulations for (50 states + DC)
v_statefips <- c('01','02','04','05','06','08','09','10','11','12','13','15','16',
                 '17','18','19','20','21','22','23','24','25','26','27','28','29',
                 '30','31','32','33','34','35','36','37','38','39','40','41','42',
                 '44','45','46','47','48','49','50','51','53','54','55','56')
