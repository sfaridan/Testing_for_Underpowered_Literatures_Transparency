
#Simulations for Testing for Underpowered Literatures
# April 10, 2025
# Stefan Faridani
# sfaridani6@gatech.edu


results_path <- paste0(root,"/output/tables/") 
figures_path <-   paste0(root,"/output/figures/") 

source(paste0(root,"/code/functions/Functions_TFUL_new_inference.R")) #Sept 2025

#Benchmark parameterization
cvs            <- c(1.96)
popsizes       <- c(500000)
cs             <- c(sqrt(2))
sigma_Ys       <- c(1)
theta0         <- c(0.9)
J_coeffs       <- c(1e-4 )
eps_coeffs     <- c(2)
noise_dgp      <- c("normal")
ns             <- c(50,500)
nsims          <- c(10000)
nu             <- c(50)
dgps           <- c("null","cauchy", "realistic", "large","worst","unif")
parms          <- expand.grid(popsizes=popsizes, cvs = cvs, cs=cs, ns=ns, sigma_Ys = sigma_Ys, theta0=theta0,J_coeffs= J_coeffs,eps_coeffs=eps_coeffs,nsims=nsims, dgps=dgps,noise_dgp=noise_dgp,nu=nu)
set.seed(1)
parms_old      <- run_sims(parms)

source(paste0(root,"/code/functions/Functions_TFUL_Feb_2026.R"))
set.seed(1)
parms_new      <- run_sims(parms)

parms_new$Cover_deltahat
parms_old$Cover_deltahat

parms_new$SD_deltahat
parms_old$SD_deltahat



