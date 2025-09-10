
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
noise_dgp      <- c("normal")
ns             <- c(50,500)
nsims          <- c(10000)
nu             <- c(50)
dgps           <- c("null","cauchy", "realistic", "large","worst","unif")
parms          <- expand.grid(popsizes=popsizes, cvs = cvs, cs=cs, ns=ns, sigma_Ys = sigma_Ys, theta0=theta0,J_coeffs= J_coeffs,eps_coeffs=eps_coeffs,nsims=nsims, dgps=dgps,noise_dgp=noise_dgp,nu=nu)
parms_out      <- run_sims(parms)
save_results(parms_out,results_path,"_Sept_25_Normal")

tab <- cbind(parms$ns,parms$dgps, 1-parms_out$beta_1, -parms_out$delta0, -parms_out$Mean_deltahat,parms_out$SD_deltahat, parms_out$SD_EST_deltahat,parms_out$Cover_deltahat)
round(tab[parms$ns==50 & noise_dgp=="normal",],2)
round(tab[parms$ns>50 & noise_dgp=="normal",],2)

noise_dgp      <- c("lognormal")
nu             <- c(185)
ns             <- c(50,500)
parms          <- expand.grid(popsizes=popsizes, cvs = cvs, cs=cs, ns=ns, sigma_Ys = sigma_Ys, theta0=theta0,J_coeffs= J_coeffs,eps_coeffs=eps_coeffs,nsims=nsims, dgps=dgps,noise_dgp=noise_dgp,nu=nu)
parms_out_non      <- run_sims(parms)
save_results(parms_out_non,results_path,"_Sept_25_lognormal")
tab <- cbind(parms_out_non$noise_dgp,parms_out_non$nu, parms_out_non$ns,parms_out_non$dgps, 1-parms_out_non$beta_1, -parms_out_non$delta0, -parms_out_non$Mean_deltahat,parms_out_non$SD_deltahat, parms_out_non$SD_EST_deltahat,parms_out_non$Cover_deltahat)
round(tab,2)


noise_dgp      <- c("t")
nu             <- c(30)
ns             <- c(50,500)
parms          <- expand.grid(popsizes=popsizes, cvs = cvs, cs=cs, ns=ns, sigma_Ys = sigma_Ys, theta0=theta0,J_coeffs= J_coeffs,eps_coeffs=eps_coeffs,nsims=nsims, dgps=dgps,noise_dgp=noise_dgp,nu=nu)
parms_out_non      <- run_sims(parms)
save_results(parms_out_non,results_path,"_Sept_25_t")
tab <- cbind(parms_out_non$noise_dgp,parms_out_non$nu, parms_out_non$ns,parms_out_non$dgps, 1-parms_out_non$beta_1, -parms_out_non$delta0, -parms_out_non$Mean_deltahat,parms_out_non$SD_deltahat, parms_out_non$SD_EST_deltahat,parms_out_non$Cover_deltahat)
round(tab,2)
