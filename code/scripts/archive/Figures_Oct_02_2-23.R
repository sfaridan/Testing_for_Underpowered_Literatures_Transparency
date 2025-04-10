#Simulations for Testing for Underpowered Literatures
# October 1st, 2023
# Stefan Faridani
# stefan.faridani@gmail.com

rm(list = ls())

library(tictoc)
library(magrittr)
library(ggplot2)
library(quadprog)
library(latex2exp)


#The normalized generalized Hermite polynomials (Carrasco 2011)
hermite_general<- function(x,j,sigma_Y){
  ans <- 0*x
  for (l in 0:floor(j/2)){
    ans <- ans+(-1)^l*factorial(2*l)/2^l/factorial(l)*choose(j,2*l)*(x/sigma_Y)^(j-2*l)
    #print((-1)^l*factorial(2*l)/2^l/factorial(l)*choose(j,2*l))
  }
  return(ans/sqrt(factorial(j)))
}

#estimate theta
estimate_theta <- function(data,cv,bandwidth){
  denominator <- mean( abs(data)<= cv+bandwidth) -  mean( abs(data)<= cv)
  numerator <- mean( abs(data)<= cv) - mean(abs(data)<= cv-bandwidth)
  return(  numerator/denominator )
}

# Estimate E[1/w_theta(T) | R]
estimate_denominator<- function(data,cv,thetahat){
  return(mean( (1+(abs(data)<cv )*(1/thetahat-1)  ) ))
}


trans_data <- function(data,J,cv,c,sigma_Y,thetahat){
  sigma_X <- sqrt(1-c^(-2)+sigma_Y^2)
  lambda1 <- sigma_Y^2/sigma_X^2
  transformed_data   <- 0*data
  for (j in seq(1,J,2)){ #only even Js will have non-zeor integ-coeff
    integrand   <- function(x){ return(hermite_general(x,j-1,sigma_X)) }
    integ_coeff <-integrate(integrand, -cv/c,cv/c)$value 
    trans_data  <- (1+(abs(data)<=cv)*(1/thetahat-1))*dnorm(data/sigma_Y)/sigma_Y*hermite_general(data,j-1,sigma_Y)
    transformed_data    <- transformed_data+ integ_coeff*(lambda1^(-(j-1)/2))*(trans_data)
  }
  return(transformed_data)
}

trans_data_delta <- function(data,J,cv,c,sigma_Y,thetahat){
  sigma_X <- sqrt(1-c^(-2)+sigma_Y^2)
  lambda1 <- sigma_Y^2/sigma_X^2
  transformed_data   <- 0*data
  for (j in seq(1,J,2)){ #only even Js will have non-zeor integ-coeff
    integrand   <- function(x){ return(hermite_general(x,j-1,sigma_X)) }
    integ_coeff <-integrate(integrand, -cv/c,cv/c)$value 
    trans_data  <- (1+(abs(data)<=cv)*(1/thetahat-1))*dnorm(data/sigma_Y)/sigma_Y*hermite_general(data,j-1,sigma_Y)
    transformed_data    <- transformed_data+ integ_coeff*(lambda1^(-(j-1)/2)-1)*(trans_data)
  }
  return(transformed_data)
}


linearize_thethat <- function(data,cv,bandwidth,thetahat){
  denom <- mean(((abs(data)<= cv)-(abs(data)<= cv-bandwidth)))
  term1 <- ((abs(data)<= cv+bandwidth)-(abs(data)<= cv))/denom
  term2 <- (1/thetahat)*((abs(data)<= cv)-(abs(data)<= cv-bandwidth))/denom
  return(term1-term2)
}

make_studymat <- function(studies){
  n <- length(studies)
  studymat <- diag(n)
  for (i in 1:n){
    for (j in i:n){
      studymat[i,j] <- studies[i] == studies[j]
      studymat[j,i] <- studymat[i,j]
    }
  }
  return(studymat)
}

get_population <- function(popsize,dgp,c){
  
  band <- 0.001
  
  #draw hs
  if(dgp=="null"){
    hs <- runif(popsize, min=-band, max=band)
  }
  else if (dgp == "Z"){
    hs <- rnorm(popsize,mean=0,sd=1)
  }
  else if (dgp == "unif"){
    hs <- runif(popsize, min=-3, max=3)
  }
  else if (dgp == "cauchy"){
    hs <- rcauchy(popsize)
  }
  else if (dgp == "on23"){
    hs <- c(rnorm(popsize,mean=2.3,sd=0.1)  )
  }
  else if (dgp == "realistic"){ #half are 80% powered, half are nearly nulls
    hs <- c(rnorm(popsize/2,mean=2.8,sd=1), rnorm(popsize/2,mean=0,sd=1)  )
  }
  else if (dgp == "large"){ #maximizes delta
    hs <- c(rnorm(popsize,mean=1.96,sd=0.2) )
  }
  else if (dgp == "worst"){ #make slope of fT at cv as negative as possible
    hs <- runif(popsize,min=1.96-1-band,max=1.96-1+band)
  }
  else{
    error(paste0(dgp," is an invalid dgp name"))
  }
  
  #convert to t-scores
  pop <- c*hs+rnorm(length(hs),mean=0,sd=1)
  return(pop)
}

estimator <- function(data,J,cv,c,sigma_Y,bandwidth,studies=NULL,include_pb=TRUE,lambda=1.45){
  output <-  list()
  n <- length(data)
  
  output$thetahat <- min(9,estimate_theta(data,cv,bandwidth))
  if(is.nan(output$thetahat) |output$thetahat<=0 ){
    output$thetahat<- 1
  }
  if(!include_pb){
    output$thetahat <- 1
  }
  gndata_num <- trans_data(data,J,cv,c,sigma_Y,output$thetahat)
  gndata_num_delta <- trans_data(data,J,cv,1,sigma_Y,output$thetahat)
  denomhat <- estimate_denominator(data,cv,output$thetahat)
  output$betahat <- mean(gndata_num)/denomhat #estimator(J,cv,c,data,sigma_Y,thetahats[sim])
  output$deltahat <- mean(gndata_num)/denomhat - mean(gndata_num_delta)/denomhat #mean(gndata_num)/denomhat-mean(abs(data)<cv)/thetahats[sim]/denomhat #mean(gndata_num_delta)/denomhat
  
  #Inference
  
  #Account for clustering within studies
  if (is.null(studies)){
    studies <- 1:n
    studymat <- diag(n)
  }
  else{
    studymat <- make_studymat(studies)
  }
  
  thlin <- linearize_thethat(data,cv,bandwidth,output$thetahat)*(1*include_pb)
  gndata <- 1*gndata_num/denomhat + thlin*mean(gndata_num*(abs(data)<cv)*(output$thetahat))/denomhat - thlin*mean((abs(data)<cv))*output$betahat/denomhat 
  gndata_delta <- gndata - (  1*gndata_num_delta/denomhat + thlin*mean(gndata_num_delta*(abs(data)<cv)*(output$thetahat))/denomhat - thlin*mean((abs(data)<cv))*( mean(gndata_num_delta)/denomhat)/denomhat ) #1*gndata_num_delta/denomhat + thlin*mean(gndata_num_delta*(abs(data)<cv)*(thetahats[sim]))/denomhat - thlin*mean((abs(data)<cv))*deltahats[sim]/denomhat 
  output$varest_theta <- ( (t(thlin)%*%studymat%*%thlin)/n  -mean(thlin)^2)/n
  output$varest_beta <- ((t(gndata)%*%studymat%*%gndata)/n -mean(gndata)^2)/n
  output$varest_delta <- ((t(gndata_delta)%*%studymat%*%gndata_delta)/n-mean(gndata_delta)^2)/n
  return(output)
}

#takes in population without publication bias and truncates it
# if absolute value of t-score is above cv, you stay
# if |T|<cv then stay with probability theta
trunc_population <- function(untrunc_population, cv, theta){
  randkeep <- runif(length(untrunc_population))
  truncated <- untrunc_population[ abs(untrunc_population)>= cv | randkeep<=theta  ]
  return(truncated)
}

run_sims<- function(parms){
  tic()
  num_parameterizations <- nrow(parms)
  
  #Outcomes
  parms$beta_c                    <- rep(NA,num_parameterizations)
  parms$beta_1                    <- rep(NA,num_parameterizations)
  parms$delta0                    <- rep(NA,num_parameterizations)
  parms$Mean_betahat              <- rep(NA,num_parameterizations)
  parms$Mean_deltahat             <- rep(NA,num_parameterizations)
  parms$Mean_thetahat             <- rep(NA,num_parameterizations)
  parms$SD_betahat                <- rep(NA,num_parameterizations)
  parms$SD_deltahat               <- rep(NA,num_parameterizations)
  parms$SD_thetahat               <- rep(NA,num_parameterizations)
  parms$SD_EST_betahat            <- rep(NA,num_parameterizations)
  parms$SD_EST_deltahat           <- rep(NA,num_parameterizations)
  parms$SD_Est_thetahat           <- rep(NA,num_parameterizations)
  parms$Cover_deltahat            <- rep(NA,num_parameterizations)
  parms$Cover_betahat             <- rep(NA,num_parameterizations)
  parms$Cover_thetahat            <- rep(NA,num_parameterizations)
  parms$J                         <- rep(NA,num_parameterizations)
  parms$epx                       <- rep(NA,num_parameterizations)
  
  
  #Loop over parameterizations
  for (parm in 1:num_parameterizations) {
    #Pre-calculate quantities for this parameterization
    set.seed(parms$seed[parm])
    
    #Draw the population of t-scores
    population_pre_pb         <- get_population(parms$popsizes[parm],parms$dgps[parm],1)
    population_counterfactual <- get_population(parms$popsizes[parm],parms$dgps[parm],parms$cs[parm])
    population_post_pb        <- trunc_population(population_pre_pb, parms$cvs[parm], parms$theta0[parm])
    
    #Record the key population estimands
    parms$beta_c[parm]        <- mean(abs(population_counterfactual) < parms$cvs[parm])
    parms$beta_1[parm]        <- mean(abs(population_pre_pb)< parms$cvs[parm])
    parms$delta0[parm]        <- parms$beta_c[parm] -  parms$beta_1[parm]
    
    #tuning parameters that depend on n
    parms$J[parm]      <- log(parms$J_coeffs[parm] * parms$ns[parm]^(-1/3)) / log(parms$sigma_Ys[parm]^2/ (1+parms$sigma_Ys[parm]^2))
    parms$eps[parm]    <- parms$eps_coeffs[parm] * (parms$ns[parm]^(-1/3))
    
    
    toc()
    tic()
    betahats       <- rep(0,nsims)
    deltahats      <- rep(0,nsims)
    thetahats      <- rep(0,nsims)
    varests        <- rep(0,nsims)
    varests_delta  <- rep(0,nsims)
    varests_theta  <- rep(0,nsims)
    for(sim in 1:parms$nsims[parm]){
      
      numsimsper <- 100
      if ( (sim) /numsimsper == round( (sim) /numsimsper) ){
        print("")
        print("")
        print(paste0("Parm: ", parm, " of ",num_parameterizations, ", Sim: ", sim, " of ", parms$nsims[parm] ))
        print(Sys.time())
        toc()
        tic()
      }
      
      #dgp
      data <- sample(x=population_post_pb,size=parms$ns[parm],replace=TRUE)
      studies <- 1:parms$ns[parm]
      
      #Estimation and inference
      output <- estimator(data,parms$J[parm],parms$cvs[parm],parms$cs[parm],parms$sigma_Ys[parm],parms$eps[parm],studies = studies)
      
      #Record results
      thetahats[sim]     <- output$thetahat
      betahats[sim]      <- output$betahat
      deltahats[sim]     <- output$deltahat
      varests[sim]       <- output$varest_beta
      varests_delta[sim] <- output$varest_delta
      varests_theta[sim] <- output$varest_theta
      print(paste0("Parm: ", parm, " of ",num_parameterizations, ", Sim: ", sim, " of ", parms$nsims[parm] ))
      #print(mean(abs(betahats[1:sim]-parms$beta1[parm])/sqrt(varests[1:sim]) <= 1.96,na.rm=1 )-mean(is.nan(varests[1:sim])))
      print(mean(abs(deltahats[1:sim]-parms$delta0[parm])/sqrt(varests_delta[1:sim]) <= 1.96,na.rm=1 )-mean(is.nan(varests_delta[1:sim])))
      #print(mean(abs(thetahats[1:sim]-parms$theta0[parm])/sqrt(varests_theta[1:sim]) <= 1.96,na.rm=1 )-mean(is.nan(varests_theta[1:sim])))
    }
    
    #Calculate and record results for this parameterization
    parms$Mean_betahat[parm]         <- mean(betahats)
    parms$Mean_deltahat[parm]        <- mean(deltahats)
    parms$Mean_thetahat[parm]        <- mean(thetahats)
    parms$SD_betahat[parm]           <- sd(betahats)
    parms$SD_deltahat[parm]          <- sd(deltahats)
    parms$SD_thetahat[parm]          <- sd(thetahats)
    parms$SD_EST_betahat[parm]       <- sqrt((mean(varests,na.rm=1)))
    parms$SD_EST_deltahat[parm]      <- sqrt((mean(varests_delta,na.rm=1)))
    parms$SD_Est_thetahat[parm]      <- sqrt((mean(varests_theta,na.rm=1)))
    parms$Cover_deltahat[parm]       <-  mean(abs(deltahats-parms$delta0[parm])/sqrt(varests_delta) <= 1.96,na.rm=1 )-mean(is.nan(varests_delta))
    parms$Cover_betahat[parm]        <-  mean(abs(betahats-parms$beta_c[parm])/sqrt(varests) <= 1.96,na.rm=1 ) -mean(is.nan(varests_delta))
    parms$Cover_thetahat[parm]       <-  mean(abs(thetahats-parms$theta0[parm])/sqrt(varests_theta) <= 1.96,na.rm=1 ) -mean(is.nan(varests_delta))
    
  }
  
  parms$RMSE_beta  <- sqrt( (parms$Mean_betahat-parms$beta_c)^2+parms$SD_deltahat^2 )
  parms$RMSE_delta <- sqrt( (parms$Mean_deltahat-parms$delta0)^2+parms$SD_betahat^2 )
  parms$RMSE_theta <- sqrt( (parms$Mean_thetahat-parms$theta0)^2+parms$SD_thetahat^2 )
  
  toc()
  return(parms)
}

#Save the results as a csv
save_results <- function(parms_out, results_path,addendum){
  sanitize_date <- gsub(as.character(Sys.Date()), pattern = "-", replacement = "_")
  sanitize_save_name <- paste0("sims_", nrow(parms_out),"_",sanitize_date)
  print( paste0(results_path, sanitize_save_name, addendum,  ".csv"))
  write.csv(x = parms_out, file = paste0(results_path, sanitize_save_name,addendum, ".csv"))
}

results_path <-  "C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/output/results/" 
figures_path <-  "C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/output/figures/" 

#Benchmark parameterization
cvs            <- c(1.96)
popsizes       <- c(100000)
cs             <- c(sqrt(2))
ns             <- c(78,500) #c(78,1000)
sigma_Ys       <- c(1)
theta0         <- c(0.5)
J_coeffs       <- c(0.001 )
eps_coeffs     <- c(1)
nsims          <- c(500)
dgps           <- c("Z","null","realistic")

#Special parameterization for this figure
eps_coeffs     <- c(0.5,1,2,8)
J_coeffs       <- sort(c(0.00005,0.0001,0.00025,0.0005,0.001,0.0025,0.005,0.05,0.025,0.01,0.1,0.2,0.3 ))
parms <- expand.grid(popsizes=popsizes, cvs = cvs, cs=cs, ns=ns, sigma_Ys = sigma_Ys, theta0=theta0,J_coeffs= J_coeffs,eps_coeffs=eps_coeffs,nsims=nsims, dgps=dgps)
parms_out <- run_sims(parms)
save_results(parms_out,results_path,"_Js_and_eps")

#Make the figures

#Read in the csv
parms_out <- read.csv("C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/output/results/sims_312_2023_10_02_Js_and_eps.csv")

parms_out <- parms_out[order(parms_out$J_coeffs),]
parms_out <- parms_out[parms_out$dgps=="null",]
e1n500 <- parms_out$eps_coeffs == 1 & parms_out$ns != 78
e05n500 <- parms_out$eps_coeffs == 0.5 & parms_out$ns != 78
e8n500 <- parms_out$eps_coeffs == 8 & parms_out$ns != 78
e2n500 <- parms_out$eps_coeffs == 2 & parms_out$ns != 78
e1n78 <- parms_out$eps_coeffs == 1 & parms_out$ns == 78
e05n78 <- parms_out$eps_coeffs == 0.5 & parms_out$ns == 78
e8n78 <- parms_out$eps_coeffs == 8 & parms_out$ns == 78
e2n78 <- parms_out$eps_coeffs == 2 & parms_out$ns == 78

xvar <- parms_out$J_coeffs 
yvar <- parms_out$Cover_deltahat

pdf(file=paste0(figures_path,"Cover_large.pdf"),width=4,height=4)
plot(xvar[e1n500],yvar[e1n500] ,ylim=c(0.4,1),log='x',xlab="D (log scale)",main="Coverage (n = 500)",ylab="Coverage", type="l")
lines(xvar[e05n500],yvar[e05n500] ,lty=c(2),col="blue")
lines(xvar[e2n500],yvar[e2n500] ,lty=c(5),col="red")
lines(xvar[e8n500],yvar[e8n500] ,lty=c(4),col="purple")
legend(x="bottomleft",legend=TeX(c("C = 0.5","C = 1", "C = 2", eps="C = 8")),col=c("blue","black","red","purple"), lty=c(2,1,5,4))
dev.off()

pdf(file=paste0(figures_path,"Cover_small.pdf"),width=4,height=4)
plot(xvar[e1n78],yvar[e1n78] ,ylim=c(0.4,1),log="x",xlab="D (log scale)",ylab="Coverage",main="Coverage (n = 78)", type="l")
lines(xvar[e05n78],yvar[e05n78] ,lty=c(2),col="blue")
lines(xvar[e2n78],yvar[e2n78] ,lty=c(5),col="red")
lines(xvar[e8n78],yvar[e8n78] ,lty=c(4),col="purple")
legend(x="bottomleft",legend=c("C = 0.5","C = 1", "C = 2", eps="C = 8"),col=c("blue","black","red","purple"), lty=c(2,1,5,4))
dev.off()

xvar <- parms_out$J_coeffs 
yvar <- parms_out$RMSE_delta
ymax <- 0.5

pdf(file=paste0(figures_path,"RMSE_large.pdf"),width=4,height=4)
plot(xvar[e1n500],yvar[e1n500] ,ylim=c(0,ymax),log="x",xlab="D (log scale)",ylab="RMSE",main="RMSE (n = 500)", type="l")
lines(xvar[e05n500],yvar[e05n500] ,lty=c(2),col="blue")
lines(xvar[e2n500],yvar[e2n500] ,lty=c(5),col="red")
lines(xvar[e8n500],yvar[e8n500] ,lty=c(4),col="purple")
legend(x="topleft",legend=c("C = 0.5","C = 1", "C = 2", eps="C = 8"),col=c("blue","black","red","purple"), lty=c(2,1,5,4))
dev.off()

pdf(file=paste0(figures_path,"RMSE_small.pdf"),width=4,height=4)
plot(xvar[e1n78],yvar[e1n78],ylim=c(0,ymax) ,log="x",xlab="D (log scale)",ylab="RMSE",main="RMSE (n = 78)", type="l")
lines(xvar[e05n78],yvar[e05n78] ,lty=c(2),col="blue")
lines(xvar[e2n78],yvar[e2n78] ,lty=c(5),col="red")
lines(xvar[e8n78],yvar[e8n78] ,lty=c(4),col="purple")
legend(x="topleft",legend=c("C = 0.5","C = 1", "C = 2", eps="C = 8"),col=c("blue","black","red","purple"), lty=c(2,1,5,4))
dev.off()

popsizes       <- c(500000)
eps_coeffs <- c(2) #c(1,2,3)
J_coeffs <- c(5e-4,1e-4) #c(1e-5,1e-4,5e-4,1e-3,5e-3)
ns <- c(50,1000)
sigma_Y <- 1
dgps <- c("null","cauchy", "realistic", "large","worst")
theta0 <- c(0.9)
parms <- expand.grid(popsizes=popsizes, cvs = cvs, cs=cs, ns=ns, sigma_Ys = sigma_Ys, theta0=theta0,J_coeffs= J_coeffs,eps_coeffs=eps_coeffs,nsims=nsims, dgps=dgps)
parms_out <- run_sims(parms)
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

#check that the empirical application still works
setwd("C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/data")
HMB <- read.csv('HummelMaedcheBenartzi.csv')
Nudge  <- readxl::read_excel("NudgeUnits.xlsx",sheet="Sheet1")

#when is this sig

Ds <- c(1e-4,2.5e-4,5e-4,7.5e-4,1e-3,2.5e-3,5e-3,7.5e-3,1e-2)
pvals <- 0*Ds
differences <- 0*Ds
ses <- 0*Ds

for (ds in 1:length(Ds)){
sigma_Y <- 1
D <- Ds[ds]
C <- 2
J_HMB <- log(D*length(length(unique(HMB$trialnumber)))^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_HMB <-  C*length(length(unique(HMB$trialnumber)))^(-1/3)
J_Nudge <- log(D*length(length(unique(Nudge$trialnumber)))^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_Nudge <-  C*length(length(unique(Nudge$trialnumber)))^(-1/3)
cv<-1.96
c <- sqrt(2)
HMB_est <- estimator(HMB$t,J_HMB,cv,c,sigma_Y,eps_HMB,studies = HMB$trialnumber)
ci_HMB <- HMB_est$deltahat +c(-1.96*sqrt(HMB_est$varest_delta),1.96*sqrt(HMB_est$varest_delta) ) 
Nudge_est <- estimator(Nudge$t,J_Nudge,cv,c,sigma_Y,eps_HMB,studies = Nudge$trialnumber)
ci_Nudge <- Nudge_est$deltahat +c(-1.96*sqrt(Nudge_est$varest_delta),1.96*sqrt(Nudge_est$varest_delta) ) 
ci_HMB
ci_Nudge
ciw_diff <- sqrt( HMB_est$varest_delta+Nudge_est$varest_delta)*1.96
diff <- HMB_est$deltahat-Nudge_est$deltahat
c(diff,ciw_diff)
pvals[ds] <- 2*(1-pnorm(abs(diff)/sqrt(HMB_est$varest_delta+Nudge_est$varest_delta)))
differences[ds] <- diff
ses[ds] <- sqrt(HMB_est$varest_delta+Nudge_est$varest_delta)
}
cbind(Ds,pvals,differences,ses)

C <- 2
D <- 5e-4
J_MM<- log(D*(length(unique(MM$article)))^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_MM <-C*(length(unique(MM$article)))^(-1/3)
MM_est <- estimator(MM$t,J_MM,cv,c,sigma_Y,eps_MM,studies = MM$article)
c(MM_est$deltahat,sqrt(MM_est$varest_delta))


