
#Simulations for Testing for Underpowered Literatures
# April 10, 2025
# Stefan Faridani
# sfaridani6@gatech.edu


table_path   <- paste0(root,"/output/tables/") 
results_path <- paste0(root,"/output/results/") 
figures_path <-   paste0(root,"/output/figures/") 

source(paste0(root,"/code/functions/Functions_TFUL.R")) 


#Benchmark normal parameterization
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
parms_normal   <- run_sims(parms)
save_results(parms_normal,results_path,"_Feb_26_Normal")
paper_tab_normal <- write_table_like_paper(parms_normal, paste0(table_path,"Table_Sim_normal_Feb26"))

####### Log normal #######
noise_dgp      <- c("lognormal")
nu             <- c(185)
ns             <- c(50,500)
parms          <- expand.grid(popsizes=popsizes, cvs = cvs, cs=cs, ns=ns, sigma_Ys = sigma_Ys, theta0=theta0,J_coeffs= J_coeffs,eps_coeffs=eps_coeffs,nsims=nsims, dgps=dgps,noise_dgp=noise_dgp,nu=nu)
results_lognorm      <- run_sims(parms)
save_results(results_lognorm,results_path,"_Feb_26_lognormal")
paper_tab_lognormal <- write_table_like_paper(results_lognorm, paste0(table_path,"Table_Sim_log_Feb26"))

######## t(30) ##########
noise_dgp      <- c("t")
nu             <- c(30)
ns             <- c(50,500)
parms          <- expand.grid(popsizes=popsizes, cvs = cvs, cs=cs, ns=ns, sigma_Ys = sigma_Ys, theta0=theta0,J_coeffs= J_coeffs,eps_coeffs=eps_coeffs,nsims=nsims, dgps=dgps,noise_dgp=noise_dgp,nu=nu)
results_t      <- run_sims(parms)
save_results(results_t,results_path,"_Feb_26_t")
paper_tab_lognormal <- write_table_like_paper(results_t, paste0(table_path,"Table_Sim_t_Feb26"))
