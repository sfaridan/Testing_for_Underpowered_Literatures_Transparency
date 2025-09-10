
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

tab <- cbind(parms_out$ns,parms_out$dgps, 1-parms_out$beta_1, -parms_out$delta0, -parms_out$Mean_deltahat,parms_out$SD_deltahat, parms_out$SD_EST_deltahat,parms_out$Cover_deltahat)
round(tab[parms_out$ns==50 & noise_dgp=="normal",],2)
round(tab[parms_out$ns>50 & noise_dgp=="normal",],2)

noise_dgp      <- c("lognormal")
nu             <- c(185)
ns             <- c(50,500)
parms          <- expand.grid(popsizes=popsizes, cvs = cvs, cs=cs, ns=ns, sigma_Ys = sigma_Ys, theta0=theta0,J_coeffs= J_coeffs,eps_coeffs=eps_coeffs,nsims=nsims, dgps=dgps,noise_dgp=noise_dgp,nu=nu)
results_lognorm      <- run_sims(parms)
save_results(results_lognorm,results_path,"_Sept_25_lognormal")
tab <- cbind(results_lognorm$noise_dgp,results_lognorm$nu, results_lognorm$ns,results_lognorm$dgps, 1-results_lognorm$beta_1, -results_lognorm$delta0, -results_lognorm$Mean_deltahat,results_lognorm$SD_deltahat, results_lognorm$SD_EST_deltahat,results_lognorm$Cover_deltahat)
saveRDS(tab, paste0(results_path,"_table_lognormal_Sept_25"))
round(tab[results_lognorm$ns==50,],2)
round(tab[results_lognorm$ns>50,],2)


noise_dgp      <- c("t")
nu             <- c(30)
ns             <- c(50,500)
parms          <- expand.grid(popsizes=popsizes, cvs = cvs, cs=cs, ns=ns, sigma_Ys = sigma_Ys, theta0=theta0,J_coeffs= J_coeffs,eps_coeffs=eps_coeffs,nsims=nsims, dgps=dgps,noise_dgp=noise_dgp,nu=nu)
results_t      <- run_sims(parms)
save_results(results_t,results_path,"_Sept_25_t")
tab <- cbind(results_t$noise_dgp,results_t$nu, results_t$ns,results_t$dgps, 1-results_t$beta_1, -results_t$delta0, -results_t$Mean_deltahat,results_t$SD_deltahat, results_t$SD_EST_deltahat,results_t$Cover_deltahat)
round(tab[results_t$ns==50,],2)
round(tab[results_t$ns>50,],2)
