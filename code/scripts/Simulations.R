
#Simulations for Testing for Underpowered Literatures
# April 10, 2025
# Stefan Faridani
# sfaridani6@gatech.edu


results_path <- paste0(root,"/output/tables/") 
figures_path <-   paste0(root,"/output/figures/") 

source(paste0(root,"/code/functions/Functions_TFUL_new_inference.R"))


#Benchmark parameterization
cvs            <- c(1.96)
popsizes       <- c(500000)
cs             <- c(sqrt(2))
sigma_Ys       <- c(1)
theta0         <- c(0.9)
J_coeffs       <- c(1e-4 )
eps_coeffs     <- c(2)
ns             <- c(50,1000)
nsims          <- c(10000)
dgps           <- c("null","cauchy", "realistic", "large","worst","unif")
parms          <- expand.grid(popsizes=popsizes, cvs = cvs, cs=cs, ns=ns, sigma_Ys = sigma_Ys, theta0=theta0,J_coeffs= J_coeffs,eps_coeffs=eps_coeffs,nsims=nsims, dgps=dgps)
parms_out      <- run_sims(parms)
save_results(parms_out,results_path,"_April_25")

tab <- cbind(parms$ns,parms$dgps, 1-parms_out$beta_1, -parms_out$delta0, -parms_out$Mean_deltahat,parms_out$SD_deltahat, parms_out$SD_EST_deltahat,parms_out$Cover_deltahat)
round(tab[parms$ns==50,],2)
round(tab[parms$ns==1000,],2)