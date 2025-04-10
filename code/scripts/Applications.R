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
cs <- c((1:20)/10 )
tc <- 1.96


#for testing
#cs <- c(1, (6:10)/5 )
#out_CR <- estimator(Tscore_camerer_rep,1:nc,tc,cs,tupper = 10,sided=2,pb=0,nsims_inf = 5)

#for testing
out_MM <- estimator(MM$t,MM$article,tc,cs,tupper = 10,sided=2,pb=1,penalty=1,nsims_inf=100)
print(cbind(cs,out_MM$ugain))
out_Nudge <- estimator(Nudge$t,Nudge$trialnumber,tc,cs,tupper = 10,sided=2,pb=1,penalty=1,nsims_inf=100)
print(cbind(cs,out_MM$ugain,out_Nudge$ugain))
out_HMB <- estimator(HMB$t,HMB$trialnumber,tc,cs,tupper = 10,sided=2,pb=1,penalty=1,nsims_inf=100)
print(cbind(cs,out_MM$ugain,out_Nudge$ugain,out_HMB$ugain))
ML_all_ts <- c(ML_Sunk$`t (equal var)`,ML_Quote$`t (equal var)`,ML_Flag$`t (equal var)`, ML_Anchoring1$`t (equal var)`, ML_Anchoring2$`t (equal var)`, ML_Anchoring3$`t (equal var)`, ML_Anchoring4$`t (equal var)`, ML_Gambler$`t (equal var)`, ML_Money$`t (equal var)`, ML_Imagined$`t (equal var)`, ML_Math$`t (equal var)` )
ML_all_sites <- c(ML_Sunk$Site,ML_Quote$Site,ML_Flag$Site, ML_Anchoring1$Site, ML_Anchoring2$Site, ML_Anchoring3$Site, ML_Anchoring4$Site, ML_Gambler$Site, ML_Money$Site, ML_Imagined$Site, ML_Math$Site )
levels(ML_all_sites) <- 1:(length(unique(ML_all_sites)))
out_ML_all <- estimator(ML_all_ts,ML_all_sites,tc,cs,tupper = 10,sided=2,pb=0,nsims_inf=100)
print(cbind(cs,out_MM$ugain,out_Nudge$ugain,out_HMB$ugain,out_ML_all$ugain))
out_CR <- estimator(Tscore_camerer_rep,1:nc,tc,cs,tupper = 10,sided=2,pb=0,nsims_inf = 100)
print(cbind(cs,out_MM$ugain,out_Nudge$ugain,out_HMB$ugain,out_ML_all$ugain,out_CR$ugain))
out_AG <- estimator(AG$t,AG$title,tc,cs,tupper = 10,sided=2,pb=1,penalty=1,nsims_inf=100)
print(cbind(cs,out_MM$ugain,out_Nudge$ugain,out_HMB$ugain,out_ML_all$ugain,out_CR$ugain,out_AG$ugain))

# #6 22 2023 1:39pm
# print(cbind(cs,out_MM$ugain,out_Nudge$ugain,out_HMB$ugain,out_ML_all$ugain,out_CR$ugain,out_AG$ugain))
# cs                           
# [1,] 0.1 0 0.00 0.00 0.00 0.00 0.00
# [2,] 0.2 0 0.00 0.00 0.47 0.00 0.00
# [3,] 0.3 0 0.00 0.00 1.00 0.00 0.00
# [4,] 0.4 0 0.00 0.01 1.00 0.00 0.00
# [5,] 0.5 0 0.04 0.01 1.00 0.00 0.00
# [6,] 0.6 0 0.17 0.01 1.00 0.01 0.00
# [7,] 0.7 0 0.26 0.01 1.00 0.03 0.00
# [8,] 0.8 0 0.33 0.02 1.00 0.06 0.00
# [9,] 0.9 0 0.43 0.02 1.00 0.11 0.00
# [10,] 1.0 0 0.00 0.00 0.00 0.00 0.00
# [11,] 1.1 1 0.46 0.98 0.00 0.77 1.00
# [12,] 1.2 1 0.42 0.98 0.00 0.73 1.00
# [13,] 1.3 1 0.38 0.98 0.00 0.69 1.00
# [14,] 1.4 1 0.35 0.98 0.00 0.65 0.99
# [15,] 1.5 1 0.32 0.98 0.00 0.60 0.99
# [16,] 1.6 1 0.26 0.98 0.00 0.56 0.98
# [17,] 1.7 1 0.25 0.98 0.00 0.49 0.98
# [18,] 1.8 1 0.23 0.98 0.00 0.44 0.98
# [19,] 1.9 1 0.21 0.98 0.00 0.40 0.97
# [20,] 2.0 1 0.20 0.97 0.00 0.36 0.97
# [21,] 2.1 1 0.18 0.96 0.00 0.31 0.96
# [22,] 2.2 1 0.17 0.95 0.00 0.30 0.94
# [23,] 2.3 1 0.16 0.95 0.00 0.26 0.92
# [24,] 2.4 1 0.16 0.95 0.00 0.24 0.91
# [25,] 2.5 1 0.15 0.94 0.00 0.23 0.89
# [26,] 2.6 1 0.15 0.94 0.00 0.23 0.87
# [27,] 2.7 1 0.13 0.93 0.00 0.22 0.83
# [28,] 2.8 1 0.12 0.90 0.00 0.19 0.81
# [29,] 2.9 1 0.12 0.89 0.00 0.17 0.77
# [30,] 3.0 1 0.11 0.88 0.00 0.16 0.77
# [31,] 3.1 1 0.11 0.86 0.00 0.16 0.75
# [32,] 3.2 1 0.11 0.82 0.00 0.15 0.74
# [33,] 3.3 1 0.11 0.82 0.00 0.14 0.71
# [34,] 3.4 1 0.11 0.81 0.00 0.13 0.71
# [35,] 3.5 1 0.11 0.77 0.00 0.12 0.66
# [36,] 3.6 1 0.11 0.76 0.00 0.12 0.64
# [37,] 3.7 1 0.11 0.76 0.00 0.12 0.62
# [38,] 3.8 1 0.10 0.73 0.00 0.11 0.59
# [39,] 3.9 1 0.10 0.72 0.00 0.11 0.59
# [40,] 4.0 1 0.10 0.71 0.00 0.09 0.58
# [41,] 4.1 1 0.10 0.69 0.00 0.09 0.58
# [42,] 4.2 1 0.10 0.68 0.00 0.09 0.56
# [43,] 4.3 1 0.10 0.68 0.00 0.09 0.55
# [44,] 4.4 1 0.10 0.67 0.00 0.09 0.54
# [45,] 4.5 1 0.10 0.66 0.00 0.09 0.51
# [46,] 4.6 1 0.09 0.65 0.00 0.08 0.51
# [47,] 4.7 1 0.08 0.64 0.00 0.08 0.51
# [48,] 4.8 1 0.07 0.64 0.00 0.08 0.50
# [49,] 4.9 1 0.07 0.62 0.00 0.08 0.48
# [50,] 5.0 1 0.07 0.61 0.00 0.08 0.46

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

makeciplot_double_gg(out_Nudge,out_HMB,  "Nudge Units", "Academic Journals",  "Nudge Intervention Experiments (DellaVigna and Linos)")
cbind(cs,out_Nudge$ubot,out_Nudge$utop)
cbind(cs,out_HMB$ubot,out_HMB$utop)

#Nudges: Academic vs Nudge Units
out_Nudge <- estimator(Nudge$t,Nudge$trialnumber,tc,cs,tupper = 10,sided=2,pb=1,penalty=1,nsims_inf=120)
out_HMB <- estimator(HMB$t,HMB$trialnumber,tc,cs,tupper = 10,sided=2,pb=1,penalty=1,nsims_inf=120)
png(filename="Nudges.png", width=900, height=500)
makeciplot_double_gg(out_Nudge,out_HMB,  "Nudge Units", "Academic Journals",  "Nudge Intervention Experiments (DellaVigna and Linos)")
dev.off()
#ggsave("output/figures/Nudges.pdf")

#Methods Matter vs AidGrade
out_MM <- estimator(MM$t,MM$article,tc,cs,tupper = 10,sided=2,pb=1,penalty=1,nsims_inf=120)
out_MM_DID <- estimator(MM_DID$t,MM_DID$article,tc,cs,tupper = 10,sided=2,pb=1,penalty=1,nsims_inf=120)
out_MM_IV <- estimator(MM_IV$t,MM_IV$article,tc,cs,tupper = 10,sided=2,pb=1,penalty=1,nsims_inf=60)
makeciplot_double_gg(out_MM,out_MM_IV,  "Brodeur RCT", "Brodeur IV",  "Brodeur")
out_AG <- estimator(AG$t,AG$title,tc,cs,tupper = 10,sided=2,pb=1,penalty=1,nsims_inf=120)
makeciplot_gg(out_MM,"Brodeur RCTs")
makeciplot_gg(out_AG,"Aidgrade RCTs")
makeciplot_double_gg(out_MM,out_AG,  "Brodeur", "AidGrade",  "Economics RCTs")
ggsave("output/figures/RCTs.pdf")

#Replications: Economics vs Psychology
out_CR <- estimator(Tscore_camerer_rep,1:nc,tc,cs,tupper = 10,sided=2,pb=0,nsims_inf = 500)
ML_all_ts <- c(ML_Sunk$`t (equal var)`,ML_Quote$`t (equal var)`,ML_Flag$`t (equal var)`, ML_Anchoring1$`t (equal var)`, ML_Anchoring2$`t (equal var)`, ML_Anchoring3$`t (equal var)`, ML_Anchoring4$`t (equal var)`, ML_Gambler$`t (equal var)`, ML_Money$`t (equal var)`, ML_Imagined$`t (equal var)`, ML_Math$`t (equal var)` )
ML_all_sites <- c(ML_Sunk$Site,ML_Quote$Site,ML_Flag$Site, ML_Anchoring1$Site, ML_Anchoring2$Site, ML_Anchoring3$Site, ML_Anchoring4$Site, ML_Gambler$Site, ML_Money$Site, ML_Imagined$Site, ML_Math$Site )
levels(ML_all_sites) <- 1:(length(unique(ML_all_sites)))
out_ML_all <- estimator(ML_all_ts,ML_all_sites,tc,cs,tupper = 10,sided=2,pb=0,nsims_inf=500)
makeciplot_double_gg(out_CR,out_ML_all,  "Economics", "Psychology",  "Systematic Laboratory Replications")
ggsave("output/figures/Replications.pdf")

#CCTs
out_CCT <- estimator(CCt_combined$t,CCt_combined$title,tc,cs,tupper = 10,sided=2,pb=1,nsims_inf = 100)
makeciplot_double_gg(out_CCT,out_MM,  "CCT", "Brodeur",  "CCTs vs all RCTs")
makeciplot_gg(out_CCT,"CCT RCTs")

out_Moral <-estimator(Moral$t,Moral$tid,tc,cs,tupper = 10,sided=2,pb=1,nsims_inf = 100)

out_McGuire <- estimator(mcGuire$t,mcGuire$title,tc,cs,tupper = 10,sided=2,pb=1,nsims_inf = 100)



out_MM_nopb <- estimator(MM$t,MM$article,tc,cs,tupper = 10,sided=2,pb=0,nsims_inf=1)



ML_MathTs <- ML_Math$`t (equal var)`[c(1:22,24:38)]
out_ML_Math <- estimator(ML_MathTs,1:length(ML_MathTs),tc,cs,tupper = 10,sided=1,pb=0)
makeciplot(out_ML_Math,main="Math Art Gender")

out_ML_Sunk <- estimator(ML_Sunk$`t (equal var)`,1:length(ML_Sunk$`t (equal var)`),tc,cs,tupper = 10,sided=1,pb=0,nsims_inf=200)
makeciplot(out_ML_Sunk,main="Sunk Costs")

out_ML_Flag <- estimator(ML_Flag$`t (equal var)`,1:length(c(ML_Flag$`t (equal var)`)),tc,cs,tupper = 10,sided=1,pb=0)
makeciplot(out_ML_Flag,main="Flag Priming")

out_CO <- estimator(Tscore_camerer_orig,1:nc,tc,cs,tupper = 10,sided=2,pb=1)
makeciplot(out_CO,main="Camerer Original")





out_AG_rcts <-  estimator(AG_rcts$t,AG_rcts$title,tc,cs,tupper = 10,sided=2,pb=1,nsims_inf=120)
makeciplot(out_AG_rcts,main="Aidgrade RCTs")
out_MM <- estimator(MM$t,MM$article,tc,cs,tupper = 10,sided=2,pb=1,nsims_inf=120)
print(out_MM$alpha_top)
makeciplot(out_MM,main="Brodeur RCTs")
makeciplot_double(out_AG_rcts,out_MM,main="RCTs",series1="AidGrade",series2="Brodeur")

out_MM_DID <- estimator(MM_DID$t,MM_DID$article,tc,cs,tupper = 10,sided=2,pb=1,nsims_inf=120)
makeciplot(out_MM_DID,main="Brodeur DIDs")
makeciplot(out_MM,main="Brodeur RCTs")
makeciplot_double(out_MM,out_MM_DID,main="Brodeur",series1="RCTs",series2="Diff in Diff")

out_MM_DD <- estimator(MM_DD$t,MM_DD$article,tc,cs,tupper = 10,sided=2,pb=1,nsims_inf=120)
makeciplot_double(out_MM,out_MM_DD,main="Brodeur",series1="RCTs",series2="RDD")
makeciplot(out_MM_DD,main="Brodeur DD")

out_MM <- estimator(MM$t,MM$article,tc,cs,tupper = 10,sided=2,pb=1,nsims_inf=120)
makeciplot(out_MM,main="Brodeur RCTs")

out_MM_pen <- estimator(MM$t,MM$article,tc,cs,tupper = 10,sided=2,pb=1,penalty=200,nsims_inf=120)
print(out_MM_pen$alpha_cibot)
makeciplot_double(out_MM,out_MM_pen,main="Brodeur",series1="no penalty",series2="penalty")

out_AG_cct <- estimator(AG_cct$t,AG_cct$title,tc,cs,tupper = 10,sided=2,pb=1)
makeciplot_double(out_AG_cct,out_MM,main="RCTs",series1="AidGrade: Deworming",series2="Brodeur")


out_AG_worm <- estimator(AG_worm$t,AG_worm$title,tc,cs,tupper = 10,sided=2,pb=1)
makeciplot(out_AG_worm,main="AG: Deworming RCTs")
makeciplot_double(out_AG_worm,out_MM,main="RCTs",series1="AidGrade",series2="Brodeur")


#report confidence intervals
print(rbind(out_MM$cs,out_MM$citop,out_MM$cibot))
print(rbind(out_CO$cs,out_CO$citop,out_CO$cibot))
print(rbind(out_CR$cs,out_CR$citop,out_CR$cibot))
print(rbind(out_AG_rcts$cs,out_AG_rcts$citop,out_AG_rcts$cibot))
