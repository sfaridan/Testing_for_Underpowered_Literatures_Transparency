rm(list = ls())

hermite_general<- function(x,j,sigma_Y){
  ans <- 0*x
  for (l in 0:floor(j/2)){
    ans <- ans+(-1)^l*factorial(2*l)/2^l/factorial(l)*choose(j,2*l)*(x/sigma_Y)^(j-2*l)
     print((-1)^l*factorial(2*l)/2^l/factorial(l)*choose(j,2*l))
    }
  return(ans/sqrt(factorial(j)))
}

estimate_theta <- function(data,cv,bandwidth){
  denominator <- mean( abs(data)<= cv+bandwidth) -  mean( abs(data)<= cv)
  numerator <- mean( abs(data)<= cv) - mean(abs(data)<= cv-bandwidth)
  return(  numerator/denominator )
}

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

get_population <- function(popsize,hgrid,hprobs){
  hprobs <- abs(hprobs)/sum(abs(hprobs))
  hs <- c()
  for (hg in 1:length(hgrid)){
    hs <- c(hs, rep(hgrid[hg],round(popsize*hprobs[hg])  ))
  }
  pop <- hs+rnorm(length(hs),mean=0,sd=1)
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
  
  output$utilitygain <- ((1+lambda)*(1-output$betahat)-1)/c^2- ((1+lambda)*(1-output$betahat+output$deltahat)-1)
  
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
  ugain_gns<-((1+lambda)*(gndata_num_delta)/denomhat)- ((1+lambda)*(gndata_num)/denomhat)/c^2
  output$varest_utilitygain <- ((t(ugain_gns)%*%studymat%*%ugain_gns)/n-mean(ugain_gns)^2)/n
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

cv <- 1.96
c <- sqrt(2)
n      <- 78
sigma_Y <- 1
theta <- 0.5
beta <- 0.1
parm <- cv-qnorm(beta)
lambda <- 1/(1-beta-(parm/2)*dnorm(cv-parm))-1
popsize <- 1000000
hgrid <- c(-0.5,1.96,1)
hprobs <- c(0,1,0)
J      <- 10*round(-log(n^(-1/4)))
bandwidth <- n^(-1/4) #try *2 for delta
nsims   <- 2000
population_pre_pb <- get_population(popsize,hgrid,hprobs)
population_counterfactual <- get_population(popsize,hgrid*c,hprobs)
beta0 <- mean(abs(population_counterfactual)< cv)
beta_pre <- mean(abs(population_pre_pb)< cv)
delta0 <- beta0 - beta_pre
population_post_pb<-trunc_population(population_pre_pb, cv, theta)
utilitystart <- ((1+lambda)*(1-beta_pre)-1)
utilityend <- ((1+lambda)*(1-beta0)-1)/2
ugain0 <- utilitystart - utilityend

#test out the estimator
betahats <- rep(0,nsims)
varests <- rep(0,nsims)
varests_delta <- rep(0,nsims)
varests_theta<- rep(0,nsims)
varests_utilitygain<- rep(0,nsims)
deltahats <- rep(0,nsims)
thetahats <- rep(0,nsims)
utilitygain <- rep(0,nsims)
numerators <- rep(0,nsims)
nones <- rep(0,nsims)
for (sim in 1:nsims){
  data <- sample(x=population_post_pb,size=n,replace=TRUE)
  
  #add some severe clustering
  studies <- c(1:length(data),1:length(data))
  data <- c(data,data)
  
  
  
  output <- estimator(data,J,cv,c,sigma_Y,bandwidth,studies = studies,lambda=lambda)
  thetahats[sim] <- output$thetahat
  betahats[sim] <- output$betahat
  deltahats[sim] <- output$deltahat
  utilitygain[sim] <- output$utilitygain
  varests[sim] <- output$varest_beta
  varests_delta[sim] <- output$varest_delta
  varests_theta[sim] <- output$varest_theta
  varests_utilitygain[sim] <- output$varest_utilitygain
  print(c("sim", sim, " of ", nsims))
  print(round(c( theta,mean(thetahats[1:sim]),sd(1/thetahats[1:sim]),sqrt(mean(varests_theta[1:sim],na.rm=1)) ),3))
  print(round(c( beta0,mean(betahats[1:sim]),sd(betahats[1:sim]) ,sqrt(mean(varests[1:sim],na.rm=1))),3))
  print(mean(abs(betahats[1:sim]-beta0)/sqrt(varests[1:sim]) > 1.96,na.rm=1 )+mean(is.nan(varests[1:sim])))
  print(round(c( delta0,mean(deltahats[1:sim],na.rm=1),sd(deltahats[1:sim],na.rm=1),sqrt(mean(varests_delta[1:sim],na.rm=1)) ),3))
  print(mean(abs(deltahats[1:sim]-delta0)/sqrt(varests_delta[1:sim]) > 1.96,na.rm=1 )+mean(is.nan(varests_delta[1:sim])))
  print(round(c( ugain0,mean(utilitygain[1:sim],na.rm=1),sd(utilitygain[1:sim],na.rm=1),sqrt(mean(varests_utilitygain[1:sim],na.rm=1)) ),3))
}
hist(betahats[1:(sim-1)])
hist(deltahats[1:(sim-1)])



#Check that empirical application works:
cv <- 1.96
c <- sqrt(2)
n      <- 78
sigma_Y <- 1
setwd("C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/data")
HMB <- read.csv('HummelMaedcheBenartzi.csv')
Nudge  <- readxl::read_excel("NudgeUnits.xlsx",sheet="Sheet1")
bandwidth_nudge <- length(Nudge$t)^(-1/4)
bandwidth_HMB <- length(HMB$t)^(-1/4)
J_nudge      <- 10*round(-log(length(Nudge$t)^(-1/4)))
J_HMB      <- 10*round(-log( length(HMB$t)^(-1/4)))
HMB_est <- estimator(HMB$t,J_HMB,cv,c,sigma_Y,bandwidth_HMB,studies = HMB$trialnumber,lambda=lambda)
ci_HMB <- HMB_est$deltahat +c(-1.96*sqrt(HMB_est$varest_delta),1.96*sqrt(HMB_est$varest_delta) ) 
Nudge_est <- estimator(Nudge$t,J_nudge,cv,c,sigma_Y,bandwidth_nudge,studies = Nudge$trialnumber,lambda=lambda)
ci_Nudge <- Nudge_est$deltahat +c(-1.96*sqrt(Nudge_est$varest_delta),1.96*sqrt(Nudge_est$varest_delta) ) 
ci_HMB
ci_Nudge
ciw_diff <- sqrt( HMB_est$varest_delta+Nudge_est$varest_delta)*1.96
diff <- HMB_est$deltahat-Nudge_est$deltahat
c(diff,ciw_diff)
c(HMB_est$thetahat,Nudge_est$thetahat)
c(HMB_est$utilitygain/sqrt(HMB_est$varest_utilitygain),Nudge_est$utilitygain/sqrt(Nudge_est$varest_utilitygain))


#Brodeur
MM = read.csv('MM_RCTs.csv')
#M<-MM[abs(MM$t)<6,] #robust to doing this
bandwidth_MM <- length(unique(MM$article))^(-1/4)*2
J_MM      <- 10*round(-log(length(unique(MM$article))^(-1/4)))
MM_est <-  estimator(MM$t,J_MM,cv,c,1,bandwidth_MM,studies = MM$article,include_pb = TRUE,lambda=lambda)
print(MM_est$utilitygain / sqrt(MM_est$varest_utilitygain))
print(MM_est$deltahat)

MM_DID = read.csv('MM_DID.csv')
bandwidth_MM_DID <- length(MM_DID$t)^(-1/4)*2
J_MM_DID      <- 10*round(-log(length(unique(MM_DID$article))^(-1/4)))
MM_DID_est <-  estimator(MM_DID$t,J_MM_DID,cv,c,1,bandwidth_MM_DID,studies = MM_DID$article,include_pb = TRUE)
print(MM_DID_est$utilitygain / sqrt(MM_DID_est$varest_utilitygain))

# pvals_rep_camerer <- c(0.16,0.012,0.001,0.003,0.571,0.0001,0.674,0.0001,0.055,0.026,0.004,0.0001,0.142,0.933,0.016,0.010,0.001,0.154)
# pvals_orig_camerer <- c(0.046,0.057,0.007,0.010, 0.033,0.001,0.01,0.001,0.03,0.011,0.001, 0.001, 0.004, 0.031,0.001,0.016,0.001,0.07)
# Tscore_camerer_orig = -qnorm(pvals_orig_camerer/2)
# Tscore_camerer_rep = -qnorm(pvals_rep_camerer/2)
# bandwidth_CO <- length(Tscore_camerer_rep)^(-1/4)*2
# J_CO     <- 10*round(-log(length(Tscore_camerer_rep)^(-1/4)))
# CO_est <-  estimator(Tscore_camerer_orig,J_CO,cv,c,sigma_Y,bandwidth_CO)
# CR_est <-  estimator(Tscore_camerer_rep,10,cv,c,1,bandwidth_CO,include_pb = TRUE)
# print(CR_est$deltahat/sqrt(CR_est$varest_delta))

#Many labs
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
ML_Math$`t (equal var)` <- as.numeric(ML_Math$`t (equal var)`)
ML_all_ts <- c(ML_Sunk$`t (equal var)`,ML_Quote$`t (equal var)`,ML_Flag$`t (equal var)`, ML_Anchoring1$`t (equal var)`, ML_Anchoring2$`t (equal var)`, ML_Anchoring3$`t (equal var)`, ML_Anchoring4$`t (equal var)`, ML_Gambler$`t (equal var)`, ML_Money$`t (equal var)`, ML_Imagined$`t (equal var)`, ML_Math$`t (equal var)` )
ML_all_sites <- c(ML_Sunk$Site,ML_Quote$Site,ML_Flag$Site, ML_Anchoring1$Site, ML_Anchoring2$Site, ML_Anchoring3$Site, ML_Anchoring4$Site, ML_Gambler$Site, ML_Money$Site, ML_Imagined$Site, ML_Math$Site )
levels(ML_all_sites) <- 1:(length(unique(ML_all_sites)))
bandwidth_ML <- length(40)^(-1/4)*2
J_ML     <- 10*round(-log((40)^(-1/4)))
out_ML_all <- estimator(ML_all_ts,J_ML,cv,c,1,bandwidth_ML,studies=ML_all_sites,include_pb = TRUE)
print((out_ML_all$deltahat)/sqrt(out_ML_all$varest_delta))
out_ML_all$deltahat+1.96*sqrt(out_ML_all$varest_delta)