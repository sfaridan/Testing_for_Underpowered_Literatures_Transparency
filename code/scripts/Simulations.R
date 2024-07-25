#Simulations for Testing for Underpowered Literatures
# July 23, 2024
# Stefan Faridani
# stefan.faridani@gmail.com

results_path <- paste0(root,"/output/results/")  #"C:/Users/stefa/OneDrive/Documents/R/Underpowered_Literatures_Transparency/output/results/" 
figures_path <- paste0(root,"/output/figures/")  #"C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/output/figures/" 

#load in functions
source(paste0(root,"/code/functions/Functions_TFUL.R"))

#Benchmark parameterization
cvs            <- c(1.96)
popsizes       <- c(500000)
cs             <- c(sqrt(2))
sigma_Ys       <- c(1)
theta0         <- c(0.9)
J_coeffs       <- c(1e-4 )
eps_coeffs     <- c(2)
ns             <- c(50,1000)
nsims          <- c(1000)
dgps           <- c("null","cauchy", "realistic", "large","worst")
parms          <- expand.grid(popsizes=popsizes, cvs = cvs, cs=cs, ns=ns, sigma_Ys = sigma_Ys, theta0=theta0,J_coeffs= J_coeffs,eps_coeffs=eps_coeffs,nsims=nsims, dgps=dgps)
parms_out      <- run_sims(parms)
save_results(parms_out,results_path,"_NEW_INFERENCE")

cbind(parms_out$Cover_deltahat, parms_out$delta0, parms_out$Mean_deltahat,parms_out$SD_deltahat, parms_out$SD_EST_deltahat)
