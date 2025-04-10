rm(list = ls())

results_path <-  "C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/output/results/" 
figures_path <-  "C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/output/figures/" 

#load in functions
source("C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/code/functions/Functions_TFUL.R")

#Benchmark parameterization
cvs            <- c(1.96)
popsizes       <- c(500000)
cs             <- c(sqrt(2))
sigma_Ys       <- c(1)
theta0         <- c(0.5)
J_coeffs       <- c(1e-4 )
eps_coeffs     <- c(1)
ns             <- c(50,1000)
nsims          <- c(2000)
dgps           <- c("null","cauchy", "realistic", "large","worst")
theta0         <- c(0.9)
parms          <- expand.grid(popsizes=popsizes, cvs = cvs, cs=cs, ns=ns, sigma_Ys = sigma_Ys, theta0=theta0,J_coeffs= J_coeffs,eps_coeffs=eps_coeffs,nsims=nsims, dgps=dgps)
parms_out      <- run_sims(parms)
save_results(parms_out,results_path,"_sims_preferred_tuning")

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