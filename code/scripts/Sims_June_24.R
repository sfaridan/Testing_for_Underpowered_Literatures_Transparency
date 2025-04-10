#Simulations for Testing for Underpowered Literatures
# October 1st, 2023
# Stefan Faridani
# stefan.faridani@gmail.com

rm(list = ls())

results_path <-  "C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/output/results/" 
figures_path <-  "C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/output/figures/" 

source("C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/code/functions/Functions_TFUL.R")

#Benchmark parameterization
cvs            <- c(1.96)
popsizes       <- c(500000)
cs             <- c(sqrt(2))
sigma_Ys       <- c(1)
theta0         <- c(0.5)
J_coeffs       <- c(1e-4,1.2e-4,.8e-4 )
eps_coeffs     <- c(2,2.4,1.8)
ns             <- c(200)
nsims          <- c(1000)
dgps           <- c("null","cauchy", "realistic", "large","worst")
parms          <- expand.grid(popsizes=popsizes, cvs = cvs, cs=cs, ns=ns, sigma_Ys = sigma_Ys, theta0=theta0,J_coeffs= J_coeffs,eps_coeffs=eps_coeffs,nsims=nsims, dgps=dgps)
parms_out      <- run_sims(parms)
#parms_out      <- run_sims_adaptive(parms)
save_results(parms_out,results_path,"_sims_preferred_tuning")

round(parms_out$Cover_deltahat,2)

#Table 1
parms_out <- read.csv("C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/output/results/sims_20_2023_10_03_sims_preferred_tuning.csv")
names(parms_out)[names(parms_out) == "Mean_betahat"] <- "Mean_bhat"
vars <- c("ns", "dgps", "beta_c","delta0", "Mean_deltahat", "SD_deltahat","SD_EST_deltahat","Cover_deltahat")
parms_out <- parms_out[order(parms_out$ns),]
picks <-  parms_out$J_coeffs == 1e-4
knitr::kable(parms_out[picks,vars], "latex",digits=2)
picks2 <-  parms_out$J_coeffs == 5e-4
knitr::kable(parms_out[picks2,vars], "latex",digits=2)
picks3 <-  parms_out$J_coeffs == 2.5e-4
knitr::kable(parms_out[picks3,vars], "latex",digits=2)

#Empirical Applications 
setwd("C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/data")
HMB <- read.csv('HummelMaedcheBenartzi.csv')
Nudge  <- readxl::read_excel("NudgeUnits.xlsx",sheet="Sheet1")
MM_DID = read.csv('MM_DID.csv')
MM = read.csv('MM_RCTs.csv')
MM_DD = read.csv('MM_DD.csv')
MM_IV = read.csv('MM_IV.csv')
MM_other = rbind(MM_DD,MM_IV,MM_DID)
AG_raw = read.csv("AidGrade-Data-v1.1 (1).csv")
AG_raw$t <- AG_raw$ES / abs(AG_raw$upper-AG_raw$lower) *1.96*2
symmetric <- abs( (AG_raw$upper-AG_raw$ES) - (AG_raw$ES-AG_raw$lower))/(abs(AG_raw$upper)+abs(AG_raw$lower)) < 0.05  
AG <- AG_raw[symmetric,]
levels(AG$title) <- 1:(length(unique(AG$title)))
titles <-  tolower(AG$title)
rcts <- grepl("rand",titles) | grepl("rct",titles) | grepl("experiment",titles) | grepl("trial",titles)
print(unique(AG$title[!rcts]))
AG_rcts <- AG[rcts,]

#Many Labs
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
ML_Math$`t (unequal var)` <- as.numeric(ML_Math$`t (unequal var)`)
ML_all_ts <- c(ML_Sunk$`t (unequal var)`,ML_Quote$`t (unequal var)`,ML_Flag$`t (unequal var)`, ML_Anchoring1$`t (unequal var)`, ML_Anchoring2$`t (unequal var)`, ML_Anchoring3$`t (unequal var)`, ML_Anchoring4$`t (unequal var)`, ML_Gambler$`t (unequal var)`, ML_Money$`t (unequal var)`, ML_Imagined$`t (unequal var)`, ML_Math$`t (unequal var)` )
ML_all_ns <- c(ML_Sunk$`t (unequal var)`,ML_Quote$`t (unequal var)`,ML_Flag$`t (unequal var)`, ML_Anchoring1$`t (unequal var)`, ML_Anchoring2$`t (unequal var)`, ML_Anchoring3$`t (unequal var)`, ML_Anchoring4$`t (unequal var)`, ML_Gambler$`t (unequal var)`, ML_Money$`t (unequal var)`, ML_Imagined$`t (unequal var)`, ML_Math$`t (unequal var)` )
ML_all_sites <- c(ML_Sunk$Site,ML_Quote$Site,ML_Flag$Site, ML_Anchoring1$Site, ML_Anchoring2$Site, ML_Anchoring3$Site, ML_Anchoring4$Site, ML_Gambler$Site, ML_Money$Site, ML_Imagined$Site, ML_Math$Site )
ML_all_ts <- ML_all_ts[ML_all_sites != "Overall:"]
ML_all_sites <- ML_all_sites[ML_all_sites != "Overall:"]
ML_all_ts <- ML_all_ts[ML_all_sites != "Mean across samples:"]
ML_all_sites <- ML_all_sites[ML_all_sites != "Mean across samples:"]
ML_all_ts <- ML_all_ts[ML_all_sites != "Overall for US participants:"]
ML_all_sites <- ML_all_sites[ML_all_sites != "Overall for US participants:"]
ML_all_ts <- ML_all_ts[ML_all_sites != "Overall (sum of samples)"]
ML_all_sites <- ML_all_sites[ML_all_sites != "Overall (sum of samples)"]

# moral behavior
moral <- read.csv("team_results competition and moral behavior.csv")
moral <- moral[moral$anap=="a",]
levels(moral$tid) <- 1:(length(unique(moral$tid)))
moral_teams <-  tolower(moral$tid)

# rpp <- read.csv("C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/data/rpp_data_CLEANED.csv")
# levels(rpp$studytitleo) <- 1:(length(unique(rpp$studytitleo)))
# rpp_titles <-  tolower(rpp$studytitleo)

Ds <- c(1e-5,1e-4)
Cs <- c(2)
cis_RCT <- estimator_wrapper(MM$t,Ds,Cs, MM$article)
cis_IV <- estimator_wrapper(MM_IV$t,Ds,Cs, MM_IV$article)
cis_DD <- estimator_wrapper(MM_DD$t,Ds,Cs, MM_DD$article)
cis_DID <- estimator_wrapper(MM_DID$t,Ds,Cs, MM_DID$article)
cis_MM_other <- estimator_wrapper(MM_other$t,Ds,Cs, MM_other$article)
cis_Nudge <- estimator_wrapper(Nudge$t,Ds,Cs,Nudge$trialnumber)
cis_ML   <- estimator_wrapper(ML_all_ts,Ds,Cs, ML_all_sites)
cis_rpp_r   <- estimator_wrapper(rpp$tstat,Ds,Cs, rpp_titles, include_pb = 0)
cis_Moral <- estimator_wrapper(moral$t,Ds,Cs, moral_teams, include_pb=0)

cis_RCT[,3:8]
cis_MM_other[,3:8]
2*(1-pnorm((cis_RCT[,3]-cis_MM_other[,3])/sqrt(cis_RCT[,4]^2+cis_MM_other[,4]^2  )))

cis_Nudge[,3:8]
cis_Moral[,3:8]
cis_ML[,3:8]


#ML exercise
(ML_Anchoring1$`Mean (High)` - ML_Anchoring1$`Mean (Low)`)/sqrt(ML_Anchoring1$`SD (High)`^2/ML_Anchoring1$`N (High)`+ML_Anchoring1$`SD (Low)`^2/ML_Anchoring1$`N (Low)`    )
ML_Anchoring1$`t (unequal var)`
