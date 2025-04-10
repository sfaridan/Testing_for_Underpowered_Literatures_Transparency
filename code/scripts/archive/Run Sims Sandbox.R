library(here)
rm(list = ls())
setwd(paste0(here(),"/R/Underpowered Literatures" ))
source("code/functions/underpowered_literatures_functions.R")
results_path <-  paste0(here(),"/R/Underpowered Literatures/output/results/" )
figures_path <-  paste0(here(),"/R/Underpowered Literatures/output/figures/" )

#define parameters
#(tgrid,hgrid,pis_true,parms$pb_param[parm],parms$tc[parm],parms$N[parm],parms$nsims[parm],parms$pb[parm],t,cs,nsims_inf=parms$nsims_inf[parm])
Mh           <- c(3)          # max h in grid
divh         <- c(10)         # one over distance between hs in grid       
Mt           <- Mh+3          # max tscore in grid
divt         <- c(120)          # one over distance between t-scores in grid
pi_dgps      <- c(6)          # which distributions pi to use (each number refers to a specific dgp)
pb_param     <- c(1)       # true probability of not publishing t-scores below t=tc
N            <- c(500)      # Number of studies in meta-sample
tc           <- c(1.96)       # t-score across which publication bias occurs
pb           <- c(0)          # indicates whether estimator should take pb into account
t            <- c(1.96)          # t-score at which we wish to evaluate counterfactual cdfs
cs           <- c(1,2) #c(0.75,0.95,1,1.05, 1.5,2,4,5.5,8)  #  counterfactual sample size increases
nsims_inf    <- c(1)        # square root of counterfactual sample size increases
nsims        <- c(100)          # Number of simulations per parameterization
sided        <- c(1)
penalty      <- c(0) #util: (1+1.45)*(diff-pwr1/2)+1/2
cv <- 1.96
beta <- 0.2
parm <- cv-qnorm(beta)
lambdas <- c(1/(1-beta-(parm/2)*dnorm(cv-parm))-1)

#run simulations 
parms <- expand.grid(lambdas=lambdas, penalty=penalty, sided=sided,nsims_inf=nsims_inf,Mh=Mh,divh=divh,Mt=Mt,divt=divt,pi_dgps=pi_dgps,pb_param=pb_param,N=N,nsims=nsims,tc=tc,pb=pb)
parms_out <- run_sims(parms,cs,t)

#record results
save_results(parms_out,results_path)

# #make figures
# [1] "Parm 1 of 5"
# [1] "Sim 20 of 100"
# [1] "Sim 40 of 100"
# [1] "Sim 60 of 100"
# [1] "Sim 80 of 100"
# [1] "Sim 100 of 100"
# [1] "cs, Increases good?, Reject incrase bad"
# cs       
# [1,] 0.75 0 0.00
# [2,] 0.95 0 0.00
# [3,] 1.00 0 0.00
# [4,] 1.05 1 1.00
# [5,] 1.50 1 1.00
# [6,] 2.00 1 1.00
# [7,] 4.00 1 0.79
# [8,] 5.50 0 0.07
# [9,] 8.00 0 0.00
# [1] "Parm 2 of 5"
# [1] "Sim 20 of 100"
# [1] "Sim 40 of 100"
# [1] "Sim 60 of 100"
# [1] "Sim 80 of 100"
# [1] "Sim 100 of 100"
# [1] "cs, Increases good?, Reject incrase bad"
# cs       
# [1,] 0.75 0 0.00
# [2,] 0.95 0 0.00
# [3,] 1.00 0 0.00
# [4,] 1.05 0 0.04
# [5,] 1.50 0 0.00
# [6,] 2.00 0 0.00
# [7,] 4.00 0 0.00
# [8,] 5.50 0 0.00
# [9,] 8.00 0 0.00
# [1] "Parm 3 of 5"
# [1] "Sim 20 of 100"
# [1] "Sim 40 of 100"
# [1] "Sim 60 of 100"
# [1] "Sim 80 of 100"
# [1] "Sim 100 of 100"
# [1] "cs, Increases good?, Reject incrase bad"
# cs    
# [1,] 0.75 0 0
# [2,] 0.95 0 0
# [3,] 1.00 0 0
# [4,] 1.05 1 1
# [5,] 1.50 1 1
# [6,] 2.00 1 1
# [7,] 4.00 1 1
# [8,] 5.50 1 1
# [9,] 8.00 1 1
# [1] "Parm 4 of 5"
# [1] "Sim 20 of 100"
# [1] "Sim 40 of 100"
# [1] "Sim 60 of 100"
# [1] "Sim 80 of 100"
# [1] "Sim 100 of 100"
# [1] "cs, Increases good?, Reject incrase bad"
# cs    
# [1,] 0.75 0 0
# [2,] 0.95 0 0
# [3,] 1.00 0 0
# [4,] 1.05 1 1
# [5,] 1.50 1 1
# [6,] 2.00 1 1
# [7,] 4.00 1 1
# [8,] 5.50 1 1
# [9,] 8.00 1 1
# [1] "Parm 5 of 5"
# [1] "Sim 20 of 100"
# [1] "Sim 40 of 100"
# [1] "Sim 60 of 100"
# [1] "Sim 80 of 100"
# [1] "Sim 100 of 100"
# [1] "cs, Increases good?, Reject incrase bad"
# cs       
# [1,] 0.75 0 0.00
# [2,] 0.95 0 0.00
# [3,] 1.00 0 0.00
# [4,] 1.05 1 1.00
# [5,] 1.50 1 1.00
# [6,] 2.00 1 1.00
# [7,] 4.00 1 1.00
# [8,] 5.50 1 0.99
# [9,] 8.00 1 0.93

#pi 6
# "cs, Increases good?, Reject incrase bad"
# cs       
# [1,] 0.01 0 0.00
# [2,] 0.10 0 0.00
# [3,] 0.25 0 0.00
# [4,] 0.50 0 0.00
# [5,] 0.75 0 0.00
# [6,] 0.80 1 0.12
# [7,] 0.99 1 1.00
# [8,] 1.00 0 0.00
# [9,] 1.01 0 0.00
# [10,] 1.50 0 0.00
