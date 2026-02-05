#Simulations for Testing for Underpowered Literatures
# October 1st, 2023
# Stefan Faridani
# stefan.faridani@gmail.com


library(tictoc)
library(magrittr)
library(ggplot2)
library(quadprog)
library(latex2exp)

de_round <- function(x){
  
  for (i in 1:length(x)){
    if(!is.na(x[i])){
      decimals <- decimalplaces(x[i])
      x[i] <- 10^(-decimals)*( 10^decimals*x[i] +(runif(1)-0.5)  ) 
    }
  }
  
  return( x )
}

#returns number of decimal places in x for de-rounding
decimalplaces <- function(x) {
  options(scipen=999)
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}


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



trans_data_new <- function(data,J,cv,c,sigma_Y){
  sigma_X <- sqrt(1-c^(-2)+sigma_Y^2)
  lambda1 <- sigma_Y^2/sigma_X^2
  transformed_data   <- 0*data
  for (j in seq(1,J,2)){ #only even Js will have non-zeor integ-coeff
    integrand   <- function(x){ return(hermite_general(x,j-1,sigma_X)) }
    integ_coeff <-(lambda1^(-(j-1)/2))*integrate(integrand, -cv/c,cv/c)$value #-integrate(integrand, -cv,cv)$value 
    trans_data  <- dnorm(data/sigma_Y)/sigma_Y*hermite_general(data,j-1,sigma_Y)
    transformed_data    <- transformed_data+ integ_coeff*(trans_data)
  }
  return(transformed_data)
}

make_studymat <- function(studies){
  studymat <- outer(studies, studies, FUN = "==")
}



make_studymat_old <- function(studies){
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

get_population <- function(popsize,dgp,c,noise_dgp="normal",nu=50){
  
  band <- 0.001
  
  #draw hs
  if(dgp=="null"){
    hs <- rep(0,popsize) #runif(popsize, min=-band, max=band)
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
  else if (dgp == "realistic"){ #half are 80% powered, half are small
    hs <- c(rnorm(popsize/2,mean=2.8,sd=1), rnorm(popsize/2,mean=0,sd=1)  )
  }
  else if (dgp == "large"){ #maximizes delta
    hs <- rnorm(popsize,mean=1.96,sd=0.2) 
  }
  else if (dgp == "worst"){ #make slope of fT at cv as negative as possible
    hs <- rnorm(popsize,mean=0.96,sd=0.2) #runif(popsize,min=1.96-1-band,max=1.96-1+band)
  }
  else{
    error(paste0(dgp," is an invalid dgp name"))
  }
  
  #convert to t-scores
  if(noise_dgp=="normal"){
    noise <- rnorm(length(hs),mean=0,sd=1)
  }else if(noise_dgp=="t"){
    noise <- rt(length(hs),nu)
  }
  else if(noise_dgp=="lognormal"){
    norms             <- rnorm(length(hs)*nu,mean=0,sd=1)
    mean_log_normal   <- exp(0.5)
    var_log_normal    <- (exp(1)-1)*exp(1)
    lognormals        <- (exp(norms)-mean_log_normal)/sqrt(var_log_normal)
    lognormals_matrix <- matrix(lognormals,nrow=length(hs),ncol=nu)
    noise             <- lognormals_matrix%*%rep( 1/sqrt(nu), nu)
  }
  pop <- c*hs+noise
  return(pop)
}

#New inference
estimator <- function(data,J,cv,c,sigma_Y,bandwidth,studies=NULL,studies2=NULL,include_pb=TRUE){
  output <-  list()
  n <- length(data)
  output$num_tscores <- length(data)
  output$num_articles <- length(unique(studies))
  
  #estimate theta
  output$thetahat <- max(min(10,estimate_theta(data,cv,bandwidth)),1/10)
  if(is.nan(output$thetahat) |output$thetahat<=0 ){
    output$thetahat<- 1
  }
  if(!include_pb){
    output$thetahat <- 1
  }
  
  
  #define useful terms
  tupper <- 1*((abs(data)>cv) * (abs(data)<= cv+bandwidth))
  tlower <- 1*((abs(data)>cv-bandwidth) * (abs(data)<= cv))
  tsmall <- 1*(abs(data)<=cv)
  Fcv <- mean( tsmall )
  thinv <- 1/output$thetahat
  
  #estimate delta

  aj_phi_sigma <- (trans_data_new(data,J,cv,c,sigma_Y)-trans_data_new(data,J,cv,1,sigma_Y))
  pb_weights <- (1+(thinv-1)*tsmall)/(1+(thinv-1)*Fcv) #removes pb
  output$deltahat <- mean(aj_phi_sigma*pb_weights)

  
  #Inference
  

  #Account for clustering within studies
  if (is.null(studies)){
    studies <- 1:n
    studymat <- diag(n)
  }
  else{
    if(is.null(studies2)){
      studymat <- make_studymat(studies)
    }
    else{
      studymat <- make_studymat(studies)+make_studymat(studies2) >0
    }
  }
  
  #variance estimation

  Xhat <- tupper/mean(tlower)-mean(tupper)/(mean(tlower)^2)*tlower
  Qhat <- mean(aj_phi_sigma*( tsmall-Fcv)/(1+Fcv*(thinv-1))^2 )
  #Bhat <-  mean(aj_phi_sigma*(thinv-1)*(1+(thinv-1)*tsmall)/(1+Fcv*(thinv-1))^2)
  Zhat <- aj_phi_sigma*(1+(thinv-1)*tsmall)/(1+(thinv-1)*Fcv)+Qhat*Xhat #-Bhat*tsmall  
  output$varest_delta <- (t(Zhat-mean(Zhat))%*%studymat%*%(Zhat-mean(Zhat)))/n^2 
  
  output$varest_theta <- (t(Xhat-mean(Xhat))%*%studymat%*%(Xhat-mean(Xhat)))/n^2 
  
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
  parms$Mean_deltahat             <- rep(NA,num_parameterizations)
  parms$Mean_thetahat             <- rep(NA,num_parameterizations)
  parms$SD_deltahat               <- rep(NA,num_parameterizations)
  parms$SD_thetahat               <- rep(NA,num_parameterizations)
  parms$SD_EST_deltahat           <- rep(NA,num_parameterizations)
  parms$SD_Est_thetahat           <- rep(NA,num_parameterizations)
  parms$Cover_deltahat            <- rep(NA,num_parameterizations)
  parms$Cover_thetahat            <- rep(NA,num_parameterizations)
  parms$J                         <- rep(NA,num_parameterizations)
  parms$epx                       <- rep(NA,num_parameterizations)
  
  
  #Loop over parameterizations
  for (parm in 1:num_parameterizations) {
    #Pre-calculate quantities for this parameterization
    set.seed(parms$seed[parm])
    
    #Draw the population of t-scores
    population_pre_pb         <- get_population(parms$popsizes[parm],parms$dgps[parm],1,noise_dgp=parms$noise_dgp[parm],nu=parms$nu[parm])
    population_counterfactual <- get_population(parms$popsizes[parm],parms$dgps[parm],parms$cs[parm],noise_dgp=parms$noise_dgp[parm],nu=parms$nu[parm])
    population_post_pb        <- trunc_population(population_pre_pb, parms$cvs[parm], parms$theta0[parm])
    
    hist(population_pre_pb)
    
    #Record the key population estimands
    parms$beta_c[parm]        <- mean(abs(population_counterfactual) < parms$cvs[parm])
    parms$beta_1[parm]        <- mean(abs(population_pre_pb)< parms$cvs[parm])
    parms$delta0[parm]        <- parms$beta_c[parm] -  parms$beta_1[parm]
    
    #tuning parameters that depend on n
    parms$J[parm]      <- log(parms$J_coeffs[parm] * parms$ns[parm]^(-1/3)) / log(parms$sigma_Ys[parm]^2/ (1+parms$sigma_Ys[parm]^2))
    parms$eps[parm]    <- parms$eps_coeffs[parm] * (parms$ns[parm]^(-1/3))
    
    
    toc()
    tic()
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
      #print(paste0("Parm: ", parm, " of ",num_parameterizations, ", Sim: ", sim, " of ", parms$nsims[parm] ))
      data <- sample(x=population_post_pb,size=parms$ns[parm],replace=TRUE)
      studies <- 1:parms$ns[parm]
  
      
      #Estimation and inference
      output <- estimator(data,parms$J[parm],parms$cvs[parm],parms$cs[parm],parms$sigma_Ys[parm],parms$eps[parm],studies = studies)
      
      #Record results
      thetahats[sim]     <- output$thetahat
      deltahats[sim]     <- output$deltahat
      varests_delta[sim] <- output$varest_delta
      
      if (sim %% 1000==0){
      
      print(paste0("Parm: ", parm, " of ",num_parameterizations, ", Sim: ", sim, " of ", parms$nsims[parm] ))
            print(mean(abs(deltahats[1:sim]-parms$delta0[parm])/sqrt(varests_delta[1:sim]) <= 1.96,na.rm=TRUE )-mean(is.nan(varests_delta[1:sim])))
      print(c((mean(varests_delta[1:sim])),sd(deltahats[1:sim])^2) )
      }
    }
    
    #Calculate and record results for this parameterization
    parms$Mean_deltahat[parm]        <- mean(deltahats)
    parms$Mean_thetahat[parm]        <- mean(thetahats)
    parms$SD_deltahat[parm]          <- sd(deltahats)
    parms$SD_thetahat[parm]          <- sd(thetahats)
    parms$SD_EST_deltahat[parm]      <- sqrt((mean(varests_delta,na.rm=TRUE)))
    parms$SD_Est_thetahat[parm]      <- sqrt((mean(varests_theta,na.rm=TRUE)))
    parms$Cover_deltahat[parm]       <-  mean(abs(deltahats-parms$delta0[parm])/sqrt(varests_delta) <= 1.96,na.rm=TRUE )-mean(is.nan(varests_delta))
    parms$Cover_thetahat[parm]       <-  mean(abs(thetahats-parms$theta0[parm])/sqrt(varests_theta) <= 1.96,na.rm=TRUE ) -mean(is.nan(varests_delta))
    
  }
  
  parms$RMSE_delta <- sqrt( (parms$Mean_deltahat-parms$delta0)^2+parms$SD_deltahat^2 )
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

estimator_wrapper <- function(data,studies,parms){
  numstudies <- length(unique(studies))
  it <- 0
  exponent <- 3
  if (include_pb==0){
    exponent <- 2
  }
  parms_out <- expand.grid(Cs= Cs, Ds = Ds, by_articles = c(0,1))
  for (parm in 1:nrows(parms_out)){
    
      
      #By no. articles
      eps <-  Cs[cc]*(numstudies)^(-1/exponent)
      J <- log(Ds[dd]*(numstudies)^(-1/exponent))/log(sigma_Y^2/(1+sigma_Y^2))
      est<- estimator(data,J,cv,c,sigma_Y,eps,studies = studies,include_pb)
      index <-  2*it-1
      cis[index,1] <- Cs[cc]
      cis[index,2] <- Ds[dd]
      cis[index,3] <-  est$deltahat 
      cis[index,4] <-  sqrt(est$varest_delta)
      cis[index,5] <-   est$deltahat-1.96*sqrt(est$varest_delta)
      cis[index,6] <-   est$deltahat+1.96*sqrt(est$varest_delta)
      cis[index,7] <-  est$thetahat 
      cis[index,8] <-  sqrt(est$varest_theta)
      
      #By no. t-scores
      eps <-  Cs[cc]*(length(data))^(-1/exponent)
      J <- log(Ds[dd]*(length(data))^(-1/exponent))/log(sigma_Y^2/(1+sigma_Y^2))
      est<- estimator(data,J,cv,c,sigma_Y,eps,studies = studies,include_pb)
      index <-  2*it
      cis[index,1] <- Cs[cc]
      cis[index,2] <- Ds[dd]
      cis[index,3] <-  est$deltahat 
      cis[index,4] <-  sqrt(est$varest_delta)
      cis[index,5] <-   est$deltahat-1.96*sqrt(est$varest_delta)
      cis[index,6] <-   est$deltahat+1.96*sqrt(est$varest_delta)
      cis[index,7] <-  est$thetahat 
      cis[index,8] <-  sqrt(est$varest_theta)
      
      print(cis[1:(index),])
  }
  return(parms_out)
}

table_values <- function(parms){
  print("Dl_hat, Std. Err, 95% CI Bot, 95% CI Top, Theta_hat, Std. Err Theta, no. Tscores, No. articles")
  print(round(c(parms$deltahat, sqrt(parms$varest_delta),parms$deltahat-1.96*sqrt(parms$varest_delta),parms$deltahat+1.96*sqrt(parms$varest_delta),parms$thetahat,sqrt(parms$varest_theta), parms$num_tscores,parms$num_articles ),3 ))
}


makeciplot_double_gg <- function(cs,deltas_1,deltas_2,ses_1,ses_2,title1, title2,main){
  level <- 1.96
  citop_1 <- abs(deltas_1)+level*ses_1
  citop_2 <- abs(deltas_2)+level*ses_2
  print(cbind(citop_1,citop_2,ses_1))
  cibot_1 <- abs(deltas_1)-level*ses_1
  cibot_2 <- abs(deltas_2)-level*ses_2
  cibot_1[cibot_1 < 0] <- 0
  cibot_2[cibot_2 < 0] <- 0
  fields <- data.frame(cs = c(cs^2,cs^2), Literature=c(rep(title1,length(cs)),rep(title2,length(cs))), citop = c(citop_1,citop_2), cibot = c(cibot_1,cibot_2),deltas = c(abs(deltas_1),abs(deltas_2)) )
  ggplot(data=fields[fields$cs >=1,], aes(x=cs, y=citop,group=Literature))  +
    theme(text = element_text(size = 20)) + 
    theme(axis.text = element_text(size = 20)) + 
    geom_ribbon(aes(ymin=cibot, ymax=citop,group=Literature,fill=Literature), linetype=2, alpha=0.4) + ggtitle(main) +
    ylab("Average Power Gain") + xlab("Sample Size Increase Factor") + 
    geom_line(aes(x=cs,y=deltas,group=Literature,color=Literature),lwd=1) + ylim(0,.3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
}

makeciplot_single_gg <- function(cs,deltas_1,ses_1,main){
  level <- 1.96
  citop_1 <- abs(deltas_1)+level*ses_1
  cibot_1 <- abs(deltas_1)-level*ses_1
  cibot_1[cibot_1 < 0] <- 0
  print(cbind(citop_1,cibot_1))
  fields <- data.frame(cs = cs^2,  citop = citop_1, cibot = cibot_1,deltas = abs(deltas_1) )
  print(fields)
  ggplot(data=fields[fields$cs >=1,], aes(x=cs, y=citop))  +
    theme(text = element_text(size = 10)) + 
    theme(axis.text = element_text(size = 10)) + 
    geom_ribbon(aes(ymin=cibot, ymax=citop), linetype=2, alpha=0.4) + ggtitle(main) +
    ylab("Average Power Gain") + xlab("Sample Size Increase Factor") + 
    geom_line(aes(x=cs,y=deltas),lwd=1) + ylim(0,.3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
}

makeciplot_triple_gg <- function(cs,deltas_1,deltas_2,deltas_3,ses_1,ses_2,ses_3,title1, title2,title3,main){
  level <- 1.96
  citop_1 <- abs(deltas_1)+level*ses_1
  citop_2 <- abs(deltas_2)+level*ses_2
  citop_3 <- abs(deltas_3)+level*ses_3
  cibot_1 <- abs(deltas_1)-level*ses_1
  cibot_2 <- abs(deltas_2)-level*ses_2
  cibot_3 <- abs(deltas_3)-level*ses_3
  cibot_1[cibot_1 < 0] <- 0
  cibot_2[cibot_2 < 0] <- 0
  cibot_3[cibot_3 < 0] <- 0
  fields <- data.frame(cs = c(cs^2,cs^2,cs^2), Literature=c(rep(title1,length(cs)),rep(title2,length(cs)),rep(title3,length(cs))), citop = c(citop_1,citop_2,citop_3), cibot = c(cibot_1,cibot_2,cibot_3),deltas = c(abs(deltas_1),abs(deltas_2),abs(deltas_3)) )
  print(fields)
  ggplot(data=fields[fields$cs >=1,], aes(x=cs, y=citop,group=Literature))  +
    theme(text = element_text(size = 10)) + 
    theme(axis.text = element_text(size = 10)) + 
    geom_ribbon(aes(ymin=cibot, ymax=citop,group=Literature,fill=Literature), linetype=2, alpha=0.4) + ggtitle(main) +
    ylab("Average Power Gain") + xlab("Sample Size Increase Factor") + 
    geom_line(aes(x=cs,y=deltas,group=Literature,color=Literature),lwd=1) + ylim(0,.3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
}



#Computes counterfactual power for a set of experiments all studying the same intervention
compute_cpower_unweighted <- function(sites, mean1, mean2, sd1, sd2, n1, n2,c=sqrt(2)){
  ind <- sites != "Overall for US participants:" & sites !="Overall:" & sites !="Mean across samples:" & sites !="Overall (sum of samples)" 
  es <- mean(mean1[ind]-mean2[ind])
  hs <- es / sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind])
  pwr_status_quo <- 1-pnorm(1.96-hs)+pnorm(-hs-1.96)
  pwr_c <- 1-pnorm(1.96-c*hs)+pnorm(-c*hs-1.96)
  rejs <- abs((mean1[ind]-mean2[ind])/sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind]))>1.96
  #uncertainty <- sqrt(sd1[1]^2/n1[1]+sd2[1]^2/n2[1]) / mean(sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind]))
  #print((mean1[ind]-mean2[ind])/sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind]))
  
  #compute standard error
  se_es <- sqrt(sum(sd1[ind]^2/n1[ind]) / length(sd1[ind])^2 +sum(sd2[ind]^2/n2[ind]) / length(sd2[ind])^2   )
  se_hs <- se_es / sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind]) 
  se_pwrs <- se_hs * (dnorm(1.96-c*hs)+dnorm(-c*hs-1.96))  #Taylor approximation
  se_pwr <- mean(se_pwrs) /sqrt(length(se_pwrs))
  print(paste0("se / es: ", se_es/es, ", se_pwr: ",se_pwr, ", delta: ",mean(pwr_status_quo)- mean(pwr_c) ))
  
  return(c(mean(rejs),mean(pwr_status_quo), mean(pwr_c), se_pwr    ))
}

compute_cpower <- function(sites, mean1, mean2, sd1, sd2, n1, n2,c=sqrt(2)){
  ind <- sites != "Overall for US participants:" & sites !="Overall:" & sites !="Mean across samples:" & sites !="Overall (sum of samples)" 
  weights <- (n1[ind]+n2[ind])/sum(n1[ind]+n2[ind])
  es <- sum((mean1[ind]-mean2[ind])*weights)
  hs <- es / sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind])
  pwr_status_quo <- 1-pnorm(1.96-hs)+pnorm(-hs-1.96)
  pwr_c <- 1-pnorm(1.96-c*hs)+pnorm(-c*hs-1.96)
  rejs <- abs((mean1[ind]-mean2[ind])/sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind]))>1.96
  #uncertainty <- sqrt(sd1[1]^2/n1[1]+sd2[1]^2/n2[1]) / mean(sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind]))
  #print((mean1[ind]-mean2[ind])/sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind]))
  
  #compute standard error
  se_es <- sqrt(sum(sd1[ind]^2/n1[ind]*weights^2) +sum(sd2[ind]^2/n2[ind]*weights^2)   )
  se_hs <- se_es / sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind]) 
  se_pwrs <- se_hs * abs(c*(dnorm(1.96-c*hs)+dnorm(-c*hs-1.96)) - (dnorm(1.96-hs)+dnorm(-hs-1.96)))  #Taylor approximation
  se_pwr <- mean(se_pwrs) /sqrt(length(se_pwrs))
  print(paste0("se / es: ", se_es/es, ", se_pwr: ",se_pwr, ", delta: ",mean(pwr_status_quo)- mean(pwr_c) ))
  
  return(c(mean(rejs),mean(pwr_status_quo), mean(pwr_c), se_pwr    ))
}

pval_comp<- function(out1,out2){
  return( 2*(1-pnorm(abs(out1$deltahat-out2$deltahat)/sqrt(out1$varest_delta+out2$varest_delta  ))))
}
