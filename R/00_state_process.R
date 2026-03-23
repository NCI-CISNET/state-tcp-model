## Load, reformat, and save state specific data as inputs for population model
## contains code for formatting and processing state specific life expectancy,
## census population data, mortality probabilities, and smoking initiation/cessation probabilities

mainDir <- "/Users/wangmengyao/Desktop/GitHub/state-tcp-model/"
setwd(file.path(mainDir))

library(reshape2)
library(readr)
library(readxl)
library(cdlTools)
library(haven)
library(dplyr)
library(stringr)
library(lubridate)
library(ggplot2)
library(ggrepel)

v_statefips=c('01','02','04','05','06','08','09',10:13,15:42,44:51,53:56)
startbc <- 1908   # starting birth cohort 
endbc <- 2100     # ending birth cohort
cohyears <- endbc-startbc+1

# Life expectancy data ----------------------------------------------------
for (fipscodeval in v_statefips){
  # Life expectancy for never smokers (NS) from 1908 to 2100
  # The files start in 1869, need to be subsetted for the correct years
  df_M.NS.LE=cbind(read.csv("data-raw/life_exp_MaleNever.csv"), 'Men', 'Never') 
  df_F.NS.LE=cbind(read.csv("data-raw/life_exp_FemaleNever.csv"), 'Women', 'Never')
  colnames(df_M.NS.LE)=colnames(df_F.NS.LE)=c('fips', 'age', 'bc', 'abbr', 'life_exp', 'sex', 'status')
  
  # Subset to years 1908–2100 only
  df_M.NS.LE = subset(df_M.NS.LE, bc >= startbc & bc <= endbc & fips == as.numeric(fipscodeval))
  df_F.NS.LE = subset(df_F.NS.LE, bc >= startbc & bc <= endbc & fips == as.numeric(fipscodeval))
  
  # Subset Life expectancy by state
  df_M.NS.LE=subset(df_M.NS.LE, fips==as.numeric(fipscodeval))
  v_M.NS.LE=df_M.NS.LE$life_exp
  m_M.NS.LE=array(v_M.NS.LE, dim=c(100,cohyears))
  
  df_F.NS.LE=subset(df_F.NS.LE, fips==as.numeric(fipscodeval))
  v_F.NS.LE=df_F.NS.LE$life_exp
  m_F.NS.LE=array(v_F.NS.LE, dim=c(100,cohyears))
  
  save(m_M.NS.LE,m_F.NS.LE,file=paste0("data/state_inputs/le_",fipscodeval,".RData"))
}


# Census population data --------------------------------------------------

for (f in v_statefips){
  df_census_data <- read_excel(paste0("data-raw/sc-est2019-syasex-",f,".xlsx"), range= "H7:AK92", col_types="numeric", col_names=FALSE)
  df_census_data= cbind(df_census_data, read_excel(paste0("data-raw/sc-est2021-syasex-",f,".xlsx"), range= "E7:J92", col_types="numeric", col_names=FALSE))
  colnames(df_census_data) = c(rep(c("B", "M", "F"), 12)) # Rename columns for male, female, both
  df_M.census_data=df_census_data[,!colnames(df_census_data) %in% c("B", "F")]
  df_F.census_data=df_census_data[,!colnames(df_census_data) %in% c("B", "M")]
  
  # Census data is from 2010-2021, so fix the population from 1908-2010 and extend forward from 2021-2100
  df_F.census_data = cbind(rep(df_F.census_data[1],102),df_F.census_data,rep(df_F.census_data[12],79))
  df_M.census_data = cbind(rep(df_M.census_data[1],102),df_M.census_data,rep(df_M.census_data[12],79))
  
  # Census data also groups all adults ages 85+ together
  # so use the SEER standard US population distribution to distribute the 85+ population by single year of age 85-99
  F.temp=df_F.census_data[86,][rep(1,15),]*c(0.163,0.143,0.126,0.106,0.091,0.077,0.064,0.053,0.042,0.034,0.028,0.021,0.015,0.011,0.026)
  df_F.census_data=rbind(df_F.census_data[1:85,], F.temp)
  colnames(df_F.census_data) = c(startbc:endbc)
  rownames(df_F.census_data)=c(0:99)
  
  M.temp=df_M.census_data[86,][rep(1,15),]*c(0.163,0.143,0.126,0.106,0.091,0.077,0.064,0.053,0.042,0.034,0.028,0.021,0.015,0.011,0.026)
  df_M.census_data=rbind(df_M.census_data[1:85,], M.temp)
  colnames(df_M.census_data) = c(startbc:endbc)
  rownames(df_M.census_data)=c(0:99)
  
  save(df_F.census_data,df_M.census_data,file=paste0("data/state_inputs/pop_",f,".RData"))
}

# Mortality data ----------------------------------------------------------
# load mortality rates 
## Current Smoker (CS), Never Smoker (NS) 2022, Former Smoker (FS) 
df_M.mortNS=read.csv('data-raw/US_mortality_state_smk_0523/all_final_results_0523_MaleNever.csv')
df_F.mortNS=read.csv('data-raw/US_mortality_state_smk_0523/all_final_results_0523_FemaleNever.csv')
df_M.mortCS=read.csv('data-raw/US_mortality_state_smk_0523/all_final_results_0523_MaleCurrent.csv')
df_F.mortCS=read.csv('data-raw/US_mortality_state_smk_0523/all_final_results_0523_FemaleCurrent.csv')
df_M.mortFS=read.csv('data-raw/US_mortality_state_smk_0523/all_final_results_0523_MaleFormer.csv')
df_F.mortFS=read.csv('data-raw/US_mortality_state_smk_0523/all_final_results_0523_FemaleFormer.csv')

df_M.mortNS$st_fips = sprintf("%02d", df_M.mortNS$st_fips)
df_F.mortNS$st_fips = sprintf("%02d", df_F.mortNS$st_fips)
df_M.mortCS$st_fips = sprintf("%02d", df_M.mortCS$st_fips)
df_F.mortCS$st_fips = sprintf("%02d", df_F.mortCS$st_fips)
df_M.mortFS$st_fips = sprintf("%02d", df_M.mortFS$st_fips)
df_F.mortFS$st_fips = sprintf("%02d", df_F.mortFS$st_fips)


for (f in v_statefips){
  
  t_init <- Sys.time() # Start timer
  cat(paste0("fips: ", f,", state:",df_M.mortNS$abbr[df_M.mortNS$st_fips==f][1]," "))  
  
  df_M.StatemortNS=subset(df_M.mortNS,st_fips==f)
  df_F.StatemortNS=subset(df_F.mortNS,st_fips==f)
  df_M.StatemortCS=subset(df_M.mortCS,st_fips==f)
  df_F.StatemortCS=subset(df_F.mortCS,st_fips==f)
  df_M.StatemortFS=subset(df_M.mortFS,st_fips==f)
  df_F.StatemortFS=subset(df_F.mortFS,st_fips==f)
  
  df_M.StatemortNS=(df_M.StatemortNS[c( "birth_year", "age","qnq","q_all","final")])
  df_F.StatemortNS=(df_F.StatemortNS[c( "birth_year", "age","qnq","q_all","final")])
  df_M.StatemortCS=(df_M.StatemortCS[c( "birth_year", "age","qnq","q_all","final")])
  df_F.StatemortCS=(df_F.StatemortCS[c( "birth_year", "age","qnq","q_all","final")])
  df_M.StatemortFS=(df_M.StatemortFS[c( "birth_year", "age","qnq","q_all","final")])
  df_F.StatemortFS=(df_F.StatemortFS[c( "birth_year", "age","qnq","q_all","final")])
  
  # extend data to 2100 keep mortality fixed after 2050
  MNSsubset=subset(df_M.StatemortNS, birth_year==max(df_M.StatemortNS$birth_year))
  for (i in (max(df_M.StatemortNS$birth_year)+1):2100){
    MNSsubset[MNSsubset==i-1]=i
    df_M.StatemortNS=rbind(df_M.StatemortNS,MNSsubset)
  }
  
  MCSsubset=subset(df_M.StatemortCS, birth_year==max(df_M.StatemortCS$birth_year))
  for (i in (max(df_M.StatemortCS$birth_year)+1):2100){
    MCSsubset[MCSsubset==i-1]=i
    df_M.StatemortCS=rbind(df_M.StatemortCS,MCSsubset)
  }
  
  MFSsubset=subset(df_M.StatemortFS, birth_year==max(df_M.StatemortFS$birth_year))
  for (i in (max(df_M.StatemortFS$birth_year)+1):2100){
    MFSsubset[MFSsubset==i-1]=i
    df_M.StatemortFS=rbind(df_M.StatemortFS,MFSsubset)
  }
  
  FNSsubset=subset(df_F.StatemortNS, birth_year==max(df_F.StatemortNS$birth_year))
  for (i in (max(df_F.StatemortNS$birth_year)+1):2100){
    FNSsubset[FNSsubset==i-1]=i
    df_F.StatemortNS=rbind(df_F.StatemortNS,FNSsubset)
  }
  
  FCSsubset=subset(df_F.StatemortCS, birth_year==max(df_F.StatemortCS$birth_year))
  for (i in (max(df_F.StatemortCS$birth_year)+1):2100){
    FCSsubset[FCSsubset==i-1]=i
    df_F.StatemortCS=rbind(df_F.StatemortCS,FCSsubset)
  }
  
  FFSsubset=subset(df_F.StatemortFS, birth_year==max(df_F.StatemortFS$birth_year))
  for (i in (max(df_F.StatemortFS$birth_year)+1):2100){
    FFSsubset[FFSsubset==i-1]=i
    df_F.StatemortFS=rbind(df_F.StatemortFS,FFSsubset)
  }
  
  # convert to matrices
  m_M.mortNS_AC=data.matrix(df_M.StatemortNS)
  m_F.mortNS_AC=data.matrix(df_F.StatemortNS)
  m_M.mortCS_AC=data.matrix(df_M.StatemortCS)
  m_F.mortCS_AC=data.matrix(df_F.StatemortCS)
  m_M.mortFS_AC=data.matrix(df_M.StatemortFS)
  m_F.mortFS_AC=data.matrix(df_F.StatemortFS)
  
  # determine mortality for FormerSmokers based on Years Since Quit
  arrlen=dim(m_M.mortNS_AC)[1]
  a_M.mortYSQ_AC=array(m_M.mortNS_AC,dim=c(arrlen,5,40))
  a_F.mortYSQ_AC=array(m_F.mortNS_AC,dim=c(arrlen,5,40))
  v_M.RRFYSQ=pmin(1,1.0313*exp(-0.024*(1:40)))
  v_F.RRFYSQ=pmin(1,0.8613*exp(-0.023*(1:40)))
  
  startbirthcoh=1908
  endyear=2100
  for (icoh in startbirthcoh:endyear){
    for (j in 1:40){
      a_M.mortYSQ_AC[100*(icoh-startbirthcoh)+(41:100),5,j]=
        pmax(v_M.RRFYSQ[j]*m_M.mortCS_AC[100*(icoh-startbirthcoh)+(41:100),5],m_M.mortNS_AC[100*(icoh-startbirthcoh)+(41:100),5])
      a_F.mortYSQ_AC[100*(icoh-startbirthcoh)+(41:100),5,j]=
        pmax(v_F.RRFYSQ[j]*m_F.mortCS_AC[100*(icoh-startbirthcoh)+(41:100),5],m_F.mortNS_AC[100*(icoh-startbirthcoh)+(41:100),5])
    }
  }
  ###### Calculate mortality by calender year #####################################
  
  colyears=2100-1908+1
  m_M.mortNS_AP=matrix(NA,nrow=100,ncol=colyears)
  m_M.mortCS_AP=matrix(NA,nrow=100,ncol=colyears)
  m_M.mortFS_AP=matrix(NA,nrow=100,ncol=colyears)
  a_M.mortYSQ_AP=array(0,dim=c(100,colyears,40))
  
  m_F.mortNS_AP=matrix(NA,nrow=100,ncol=colyears)
  m_F.mortCS_AP=matrix(NA,nrow=100,ncol=colyears)
  m_F.mortFS_AP=matrix(NA,nrow=100,ncol=colyears)
  a_F.mortYSQ_AP=array(0,dim=c(100,colyears,40))
  
  ### we only have data as early as 1908 so fix data before then
  calstart=1908
  datastart=1908
  for (i in calstart:endyear){
    for (age in 0:99){
      byr=i-age
      if(byr<calstart){
        m_M.mortNS_AP[age+1,i-(calstart-1)]=m_M.mortNS_AC[(m_M.mortNS_AC[,1]==datastart)&(m_M.mortNS_AC[,2]==age),5]
        m_M.mortCS_AP[age+1,i-(calstart-1)]=m_M.mortCS_AC[(m_M.mortCS_AC[,1]==datastart)&(m_M.mortCS_AC[,2]==age),5]
        m_M.mortFS_AP[age+1,i-(calstart-1)]=m_M.mortFS_AC[(m_M.mortFS_AC[,1]==datastart)&(m_M.mortFS_AC[,2]==age),5]
        m_F.mortNS_AP[age+1,i-(calstart-1)]=m_F.mortNS_AC[(m_F.mortNS_AC[,1]==datastart)&(m_F.mortNS_AC[,2]==age),5]
        m_F.mortCS_AP[age+1,i-(calstart-1)]=m_F.mortCS_AC[(m_F.mortCS_AC[,1]==datastart)&(m_F.mortCS_AC[,2]==age),5]
        m_F.mortFS_AP[age+1,i-(calstart-1)]=m_F.mortFS_AC[(m_F.mortFS_AC[,1]==datastart)&(m_F.mortFS_AC[,2]==age),5]
      }
      else{
        m_M.mortNS_AP[age+1,i-(calstart-1)]=m_M.mortNS_AC[(m_M.mortNS_AC[,1]==byr)&(m_M.mortNS_AC[,2]==age),5]
        m_M.mortCS_AP[age+1,i-(calstart-1)]=m_M.mortCS_AC[(m_M.mortCS_AC[,1]==byr)&(m_M.mortCS_AC[,2]==age),5]
        m_M.mortFS_AP[age+1,i-(calstart-1)]=m_M.mortFS_AC[(m_M.mortFS_AC[,1]==byr)&(m_M.mortFS_AC[,2]==age),5]
        m_F.mortNS_AP[age+1,i-(calstart-1)]=m_F.mortNS_AC[(m_F.mortNS_AC[,1]==byr)&(m_F.mortNS_AC[,2]==age),5]
        m_F.mortCS_AP[age+1,i-(calstart-1)]=m_F.mortCS_AC[(m_F.mortCS_AC[,1]==byr)&(m_F.mortCS_AC[,2]==age),5]
        m_F.mortFS_AP[age+1,i-(calstart-1)]=m_F.mortFS_AC[(m_F.mortFS_AC[,1]==byr)&(m_F.mortFS_AC[,2]==age),5]
      }
      
      for (j in 1:40){
        if(byr<calstart){
          a_M.mortYSQ_AP[age+1,i-(calstart-1),j]=a_M.mortYSQ_AC[(a_M.mortYSQ_AC[,1,j]==datastart)&(a_M.mortYSQ_AC[,2,j]==age),5,j]
          a_F.mortYSQ_AP[age+1,i-(calstart-1),j]=a_F.mortYSQ_AC[(a_F.mortYSQ_AC[,1,j]==datastart)&(a_F.mortYSQ_AC[,2,j]==age),5,j]}
        else{
          a_M.mortYSQ_AP[age+1,i-(calstart-1),j]=a_M.mortYSQ_AC[(a_M.mortYSQ_AC[,1,j]==byr)&(a_M.mortYSQ_AC[,2,j]==age),5,j]
          a_F.mortYSQ_AP[age+1,i-(calstart-1),j]=a_F.mortYSQ_AC[(a_F.mortYSQ_AC[,1,j]==byr)&(a_F.mortYSQ_AC[,2,j]==age),5,j]
        }
      }
    }
  }
  
  save(m_M.mortNS_AC,m_F.mortNS_AC,
       m_M.mortCS_AC,m_F.mortCS_AC,
       m_M.mortFS_AC,m_F.mortFS_AC, 
       a_M.mortYSQ_AC, a_F.mortYSQ_AC,
       m_M.mortNS_AP,m_M.mortCS_AP,m_M.mortFS_AP,a_M.mortYSQ_AP,
       m_F.mortNS_AP,m_F.mortCS_AP,m_F.mortFS_AP,a_F.mortYSQ_AP,
       file=paste0("data/state_inputs/mort_rates/mort_",f,".RData"))
  
  print(Sys.time() - t_init) # End timer
}


# Load in mortality rates and convert them into probabilities 

for (f in v_statefips){
  
  load(paste0("data/state_inputs/mort_rates/mort_",f,".RData")) #mortality
  
  # Define the function to compute the probability
  compute_probability <- function(rate) {
    t <- 1  # time in years
    p <- 1 - exp(-rate * t)
    return(p)
  }
  
  m_p_M.mortNS_AC <- m_M.mortNS_AC
  m_p_F.mortNS_AC <- m_F.mortNS_AC
  m_p_M.mortCS_AC <- m_M.mortCS_AC
  m_p_F.mortCS_AC <- m_F.mortCS_AC
  m_p_M.mortFS_AC <- m_M.mortFS_AC
  m_p_F.mortFS_AC <- m_F.mortFS_AC
  a_p_M.mortYSQ_AC <- a_M.mortYSQ_AC
  a_p_F.mortYSQ_AC <- a_F.mortYSQ_AC
  
  m_p_M.mortNS_AC[, 5] <- compute_probability(m_M.mortNS_AC[, 5])
  m_p_F.mortNS_AC[, 5] <- compute_probability(m_F.mortNS_AC[, 5])
  m_p_M.mortCS_AC[, 5] <- compute_probability(m_M.mortCS_AC[, 5])
  m_p_F.mortCS_AC[, 5] <- compute_probability(m_F.mortCS_AC[, 5])
  m_p_M.mortFS_AC[, 5] <- compute_probability(m_M.mortFS_AC[, 5])
  m_p_F.mortFS_AC[, 5] <- compute_probability(m_F.mortFS_AC[, 5])
  a_p_M.mortYSQ_AC[, 5,] <- compute_probability(a_M.mortYSQ_AC[, 5,])
  a_p_F.mortYSQ_AC[, 5,] <- compute_probability(a_F.mortYSQ_AC[, 5,])
  
  # Apply the function to all columns for AP files
  m_p_M.mortNS_AP <- compute_probability(m_M.mortNS_AP)
  m_p_M.mortCS_AP <- compute_probability(m_M.mortCS_AP)
  m_p_M.mortFS_AP <- compute_probability(m_M.mortFS_AP)
  m_p_F.mortNS_AP <- compute_probability(m_F.mortNS_AP)
  m_p_F.mortCS_AP <- compute_probability(m_F.mortCS_AP)
  m_p_F.mortFS_AP <- compute_probability(m_F.mortFS_AP)
  a_p_F.mortYSQ_AP <- apply(a_F.mortYSQ_AP, c(1, 2, 3), compute_probability)
  a_p_M.mortYSQ_AP <- apply(a_M.mortYSQ_AP, c(1, 2, 3), compute_probability)
  
  
  save(
    m_p_M.mortNS_AC, m_p_F.mortNS_AC,
    m_p_M.mortCS_AC, m_p_F.mortCS_AC,
    m_p_M.mortFS_AC, m_p_F.mortFS_AC,
    a_p_M.mortYSQ_AC, a_p_F.mortYSQ_AC,
    m_p_M.mortNS_AP, m_p_M.mortCS_AP, m_p_M.mortFS_AP, a_p_M.mortYSQ_AP,
    m_p_F.mortNS_AP, m_p_F.mortCS_AP, m_p_F.mortFS_AP, a_p_F.mortYSQ_AP,
    file = paste0("data/state_inputs/mort_rates/p.mort_", f, ".RData"))
  
}


# Smoking inputs ----------------------------------------------------------
# state initiation and cessation probabilities, load and reformat by cohort

# Read in smoking parameters for all states provided by Ted
df_M.params=read.csv('data-raw/params_022422_1.csv') 
df_F.params=read.csv('data-raw/params_022422_2.csv')

for (f in v_statefips){
  
  t_init <- Sys.time() # Start timer
  
  df_M.state =subset(df_M.params,st_fips==as.numeric(f))
  df_F.state =subset(df_F.params,st_fips==as.numeric(f))
  
  df_M.stateinit=df_M.state[c("age","coh","per","init")]
  df_M.statecess=df_M.state[c("age","coh","per","quit")]
  df_M.statenever=df_M.state[c("age","coh","per","nver")]
  df_M.stateformer=df_M.state[c("age","coh","per","form")]
  df_M.statecurrent=df_M.state[c("age","coh","per","curr")]
  
  df_F.stateinit=df_F.state[c("age","coh","per","init")]
  df_F.statecess=df_F.state[c("age","coh","per","quit")]
  df_F.statenever=df_F.state[c("age","coh","per","nver")]
  df_F.statecurrent=df_F.state[c("age","coh","per","curr")]
  df_F.stateformer=df_F.state[c("age","coh","per","form")]
  
  df_M.initmelt=melt(df_M.stateinit,id=c("age","coh","per"))
  df_M.cessmelt=melt(df_M.statecess,id=c("age","coh","per"))
  df_M.nevermelt=melt(df_M.statenever,id=c("age","coh","per"))
  df_M.currentmelt=melt(df_M.statecurrent,id=c("age","coh","per"))
  df_M.formermelt=melt(df_M.stateformer,id=c("age","coh","per"))
  
  df_F.initmelt=melt(df_F.stateinit,id=c("age","coh","per"))
  df_F.cessmelt=melt(df_F.statecess,id=c("age","coh","per"))
  df_F.nevermelt=melt(df_F.statenever,id=c("age","coh","per"))
  df_F.currentmelt=melt(df_F.statecurrent,id=c("age","coh","per"))
  df_F.formermelt=melt(df_F.stateformer,id=c("age","coh","per"))
  
  # AC 'age-cohort'
  df_M.initAC=dcast(df_M.initmelt,age~coh)
  df_M.cessAC=dcast(df_M.cessmelt,age~coh)
  df_M.neverAC=dcast(df_M.nevermelt,age~coh)
  df_M.currentAC=dcast(df_M.currentmelt,age~coh)
  df_M.formerAC=dcast(df_M.formermelt,age~coh)
  
  # convert to matrix
  m_M.initAC = as.matrix(df_M.initAC[-c(1)])
  m_M.cessAC = as.matrix(df_M.cessAC[-c(1)])
  m_M.neverAC = as.matrix(df_M.neverAC[-c(1)])
  m_M.currentAC = as.matrix(df_M.currentAC[-c(1)])
  m_M.formerAC = as.matrix(df_M.formerAC[-c(1)])
  
  # Female to 'age-cohort'
  df_F.initAC=dcast(df_F.initmelt,age~coh)
  df_F.cessAC=dcast(df_F.cessmelt,age~coh)
  df_F.neverAC=dcast(df_F.nevermelt,age~coh)
  df_F.currentAC=dcast(df_F.currentmelt,age~coh)
  df_F.formerAC=dcast(df_F.formermelt,age~coh)
  
  # convert to matrix
  m_F.initAC = as.matrix(df_F.initAC[-c(1)])
  m_F.cessAC = as.matrix(df_F.cessAC[-c(1)])
  m_F.neverAC = as.matrix(df_F.neverAC[-c(1)])
  m_F.currentAC = as.matrix(df_F.currentAC[-c(1)])
  m_F.formerAC = as.matrix(df_F.formerAC[-c(1)])
  
  save(m_M.initAC, m_M.cessAC, m_M.neverAC, m_M.currentAC, m_M.formerAC,
       m_F.initAC, m_F.cessAC, m_F.neverAC, m_F.currentAC, m_F.formerAC, 
       file=paste0("data/state_inputs/smk_",f,".RData"))
  print(Sys.time() - t_init) # End timer
}

## Load data

for (f in v_statefips) {
  load(paste0("data/state_inputs/le_",f,".RData"))
  load(paste0("data/state_inputs/pop_",f,".RData"))
  load(paste0("data/state_inputs/mort_rates/mort_",f,".RData"))
  load(paste0("data/state_inputs/mort_rates/p.mort_", f, ".RData"))
  load(paste0("data/state_inputs/smk_",f,".RData"))
}


## Sources:
# Surveillance Epidemiology and End Results (SEER) Program. Standard Populations - Single Ages. Accessed 2/19/2024, 2024. https://seer.cancer.gov/stdpopulations/stdpop.singleages.html
