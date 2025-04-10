rm(list = ls())
library(here)
library(plotrix) 
library(ggplot2)
setwd(paste0(here(),"/R/Underpowered Literatures" ))
source("code/functions/underpowered_literatures_functions.R")

#Many Labs
setwd(paste0(here(),"/R/Underpowered Literatures/data" ))
ML_Anchoring1 = readxl::read_excel("MLSS.xlsx",sheet="Anchoring1")
ML_Anchoring2 = readxl::read_excel("MLSS.xlsx",sheet="Anchoring2")
ML_Anchoring3 = readxl::read_excel("MLSS.xlsx",sheet="Anchoring3")
ML_Anchoring4 = readxl::read_excel("MLSS.xlsx",sheet="Anchoring4")
ML_Gambler = readxl::read_excel("MLSS.xlsx",sheet="Gambler's Fallacy")
ML_Money = readxl::read_excel("MLSS.xlsx",sheet="Money Priming")
ML_Imagined = readxl::read_excel("MLSS.xlsx",sheet="Imagined Contact")
ML_Sunk = readxl::read_excel("MLSS.xlsx",sheet="Sunk Costs")
ML_Quote = readxl::read_excel("MLSS.xlsx",sheet="Quote Attribution")
ML_Flag = readxl::read_excel("MLSS.xlsx",sheet="Flag Priming")
ML_Math = readxl::read_excel("MLSS.xlsx",sheet="Math_Art Gender")
ML_Math <- ML_Math[c(1:22,24:28),]
ML_Math$`t (equal var)` <- as.numeric(ML_Math$`t (equal var)`)

#Moral behavior
Moral = read.csv(paste0(here(),'/R/Underpowered Literatures/data/team_results competition and moral behavior.csv'))
levels(Moral$tid) <- 1:(length(unique(Moral$tid)))


#DellaVigna and Linos
#HMB
HMB = read.csv(paste0(here(),'/R/Underpowered Literatures/data/HummelMaedcheBenartzi.csv'))
Nudge  = readxl::read_excel("NudgeUnits.xlsx",sheet="Sheet1")



#Methods Matter
MM = read.csv(paste0(here(),'/R/Underpowered Literatures/data/MM_RCTs.csv'))
MM_DID = read.csv(paste0(here(),'/R/Underpowered Literatures/data/MM_DID.csv'))
MM_DD = read.csv(paste0(here(),'/R/Underpowered Literatures/data/MM_DD.csv'))
MM_IV = read.csv(paste0(here(),'/R/Underpowered Literatures/data/MM_IV.csv'))


#Camerer (copied from table in paper)
pvals_rep_camerer <- c(0.16,0.012,0.001,0.003,0.571,0.0001,0.674,0.0001,0.055,0.026,0.004,0.0001,0.142,0.933,0.016,0.010,0.001,0.154)
pvals_orig_camerer <- c(0.046,0.057,0.007,0.010, 0.033,0.001,0.01,0.001,0.03,0.011,0.001, 0.001, 0.004, 0.031,0.001,0.016,0.001,0.07)
Tscore_camerer_orig = -qnorm(pvals_orig_camerer/2)
Tscore_camerer_rep = -qnorm(pvals_rep_camerer/2)

#negative signs for replications
Tscore_camerer_rep[7] = -Tscore_camerer_rep[7]
Tscore_camerer_rep[14] = -Tscore_camerer_rep[14]
Tscore_camerer_rep[18] = -Tscore_camerer_rep[18]



nc <- length(Tscore_camerer_orig)

#AidGrade
#get Aidgrade data
setwd(paste0(here(),"/R/Underpowered Literatures" ))
AG_raw = read.csv("data/AidGrade-Data-v1.1 (1).csv")

#compute T scores from CIs
AG_raw$t <- AG_raw$ES / abs(AG_raw$upper-AG_raw$lower) *1.96*2

#delete 20 rows with ci's that don't make sense
symmetric <- abs( (AG_raw$upper-AG_raw$ES) - (AG_raw$ES-AG_raw$lower))/(abs(AG_raw$upper)+abs(AG_raw$lower)) < 0.05  
AG <- AG_raw[symmetric,]

#give articles numbers
levels(AG$title) <- 1:(length(unique(AG$title)))

#keep only rcts
titles <-  tolower(AG$title)
rcts <- grepl("rand",titles) | grepl("rct",titles) | grepl("experiment",titles) | grepl("trial",titles)
print(unique(AG$title[!rcts]))
AG_rcts <- AG[rcts,]
AG_cct <- AG_rcts[ AG_rcts$intervention=="Conditional Cash Transfers",]
AG_worm <- AG_rcts[ AG_rcts$intervention=="Deworming",]


#Mcguire CCTs
mcGuire = read.csv("data/McGuire CCT MetaStudy.csv")
CCt_combined = data.frame( t= c(AG_cct$t,mcGuire$t ),title= c(AG_cct$title,mcGuire$title ) )
levels(CCt_combined$title) <- 1:(length(unique(CCt_combined$title)))
levels(mcGuire$title) <- 1:(length(unique(mcGuire$title)))


#Define research questions
cs <- c((1:50)/10 )
tc <- 1.96


#for testing
#cs <- c(1, (6:10)/5 )
#out_CR <- estimator(Tscore_camerer_rep,1:nc,tc,cs,tupper = 10,sided=2,pb=0,nsims_inf = 5)

#for testing
out_MM <- estimator(MM$t,MM$article,tc,cs,tupper = 10,sided=2,pb=1,penalty=1,nsims_inf=100)
print(c(out_MM$ubot,out_MM$utop))
out_Nudge <- estimator(Nudge$t,Nudge$trialnumber,tc,cs,tupper = 10,sided=2,pb=1,penalty=1,nsims_inf=100)
print(c(out_Nudge$ubot,out_Nudge$utop))
out_HMB <- estimator(HMB$t,HMB$trialnumber,tc,cs,tupper = 10,sided=2,pb=1,penalty=1,nsims_inf=100)
print(c(out_MM$ubot,out_MM$utop))
print(c(out_Nudge$ubot,out_Nudge$utop))
print(c(out_HMB$ubot,out_HMB$utop))
# 6 21 23 3:01pm
#print(c(out_MM$ubot,out_MM$utop))
#2.5% 97.5% 
#  3     5 
#> print(c(out_Nudge$ubot,out_Nudge$utop))
#2.5% 97.5% 
#  0.7   3.0 
#> print(c(out_HMB$ubot,out_HMB$utop))
#2.5%  97.5% 
#  1.3475 3.4525 


#print(cbind(parms_out$N, parms_out$cover_optim,parms_out$optim,parms_out$rmse_optim) )
#[,1] [,2] [,3]      [,4]
#[1,]   18 0.86  2.1 1.6704790
#[2,]   18 0.96  3.9 0.9519979
#[3,]   80 0.93  2.1 1.7445343
#[4,]   80 0.95  3.9 0.7379024
#[5,] 2000 0.80  2.1 0.7577599
#[6,] 2000 0.93  3.9 0.6926760

#622
#,parms_out$rmse_optim) )
#[,1] [,2] [,3]      [,4]
#[1,]   18 0.77  2.1 1.8781906
#[2,]   18 0.95  3.9 0.8785784
#[3,]   80 0.84  2.1 1.7024982
#[4,]   80 0.89  3.9 0.7395269
#[5,] 2000 0.96  2.1 0.6442049
#[6,] 2000 0.96  3.9 0.6146544

[1,] 0.50 0 0.00
[2,] 0.75 0 0.00
[3,] 1.00 0 0.00
[4,] 1.50 1 0.99
[5,] 2.00 1 0.91
[6,] 4.00 1 0.19
[7,] 5.50 0 0.02
     8.00 0 0.00
> 
