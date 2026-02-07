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



deconvolve_data <- function(data,J,cv,c,sigma_Y){
  sigma_X <- sqrt(1-c^(-2)+sigma_Y^2)
  lambda1 <- sigma_Y^2/sigma_X^2
  transformed_data   <- 0*data
  for (j in seq(0,J,2)){ #only even js will have non-zero aj
    phi_j               <- function(x){ return(hermite_general(x,j,sigma_X)) }
    psi_j               <- function(x){ return(hermite_general(x,j,sigma_Y)) }
    aj                  <- integrate(psi_j, -cv,cv)$value- (lambda1^(-(j)/2))*integrate(phi_j, -cv/c,cv/c)$value 
    trans_data_j        <- aj*dnorm(data/sigma_Y)/sigma_Y*hermite_general(data,j,sigma_Y)
    transformed_data    <- transformed_data + trans_data_j
  }
  return(transformed_data)
}


make_studymat <- function(studies){
  studymat <- outer(studies, studies, FUN = "==")
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

  data_deconvolved  <- deconvolve_data(data,J,cv,c,sigma_Y) #-deconvolve_data(data,J,cv,1,sigma_Y)) #redo this to be clearer
  pb_weights        <- (1+(thinv-1)*tsmall)/(1+(thinv-1)*Fcv) #removes pb 
  output$deltahat   <- mean(data_deconvolved*pb_weights)

  
  #Inference
  

  #Account for clustering within studies
  # studies2: for two-way clustering
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
  Qhat <- mean(data_deconvolved*( tsmall-Fcv)/(1+Fcv*(thinv-1))^2 )
  Zhat <- data_deconvolved*(1+(thinv-1)*tsmall)/(1+(thinv-1)*Fcv)+Qhat*Xhat 
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
    parms$delta0[parm]        <- parms$beta_1[parm]- parms$beta_c[parm] 
    
    #tuning parameters that depend on n
    parms$J[parm]      <- log(parms$J_coeffs[parm] * parms$ns[parm]^(-1/3)) / log(parms$sigma_Ys[parm]^2/ (1+parms$sigma_Ys[parm]^2))
    parms$eps[parm]    <- parms$eps_coeffs[parm] * (parms$ns[parm]^(-1/3))
    
    
    toc()
    tic()
    deltahats      <- rep(0,parms$nsims[parm])
    thetahats      <- rep(0,parms$nsims[parm])
    varests        <- rep(0,parms$nsims[parm])
    varests_delta  <- rep(0,parms$nsims[parm])
    varests_theta  <- rep(0,parms$nsims[parm])
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
      varests_theta[sim] <- output$varest_theta
      
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

table_values <- function(parms){
  print("Dl_hat, Std. Err, 95% CI Bot, 95% CI Top, Theta_hat, Std. Err Theta, no. Tscores, No. articles")
  print(round(c(parms$deltahat, sqrt(parms$varest_delta),parms$deltahat-1.96*sqrt(parms$varest_delta),parms$deltahat+1.96*sqrt(parms$varest_delta),parms$thetahat,sqrt(parms$varest_theta), parms$num_tscores,parms$num_articles ),3 ))
}


makeciplot_double_gg <- function(cs,deltas_1,deltas_2,ses_1,ses_2,title1, title2,main){
  level <- 1.96
  citop_1 <- (deltas_1)+level*ses_1
  citop_2 <- (deltas_2)+level*ses_2
  print(cbind(citop_1,citop_2,ses_1))
  cibot_1 <- (deltas_1)-level*ses_1
  cibot_2 <- (deltas_2)-level*ses_2
  cibot_1[cibot_1 < 0] <- 0
  cibot_2[cibot_2 < 0] <- 0
  fields <- data.frame(cs = c(cs^2,cs^2), Literature=c(rep(title1,length(cs)),rep(title2,length(cs))), citop = c(citop_1,citop_2), cibot = c(cibot_1,cibot_2),deltas = c((deltas_1),(deltas_2)) )
  ggplot(data=fields[fields$cs >=1,], aes(x=cs, y=citop,group=Literature))  +
    theme(text = element_text(size = 20)) + 
    theme(axis.text = element_text(size = 20)) + 
    geom_ribbon(aes(ymin=cibot, ymax=citop,group=Literature,fill=Literature), linetype=2, alpha=0.4) + ggtitle(main) +
    ylab("Average Power Gain") + xlab("Sample Size Increase Factor") + 
    geom_line(aes(x=cs,y=deltas,group=Literature,color=Literature),lwd=1) + ylim(0,.3) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
}

makeciplot_single_gg <- function(cs,deltas_1,ses_1,main){
  level <- 1.96
  citop_1 <- (deltas_1)+level*ses_1
  cibot_1 <- (deltas_1)-level*ses_1
  cibot_1[cibot_1 < 0] <- 0
  print(cbind(citop_1,cibot_1))
  fields <- data.frame(cs = cs^2,  citop = citop_1, cibot = cibot_1,deltas = (deltas_1) )
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
  fields <- data.frame(cs = c(cs^2,cs^2,cs^2), Literature=c(rep(title1,length(cs)),rep(title2,length(cs)),rep(title3,length(cs))), citop = c(citop_1,citop_2,citop_3), cibot = c(cibot_1,cibot_2,cibot_3),deltas = c((deltas_1),(deltas_2),(deltas_3)) )
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



write_table_like_paper <- function(df, tex_file,
                                   dgp_order = c("null","cauchy","realistic","large","worst","unif")) {
  
  # Map your internal DGP names to the paper's labels
  dgp_label <- function(x) {
    # IMPORTANT: adjust these labels ONLY if your paper uses different names.
    out <- x
    out[out == "null"]      <- "True Null"
    out[out == "cauchy"]    <- "Cauchy"
    out[out == "realistic"] <- "Bimodal"
    out[out == "large"]     <- "Large"
    out[out == "worst"]     <- "Slope"
    out[out == "unif"]      <- "Uniform"
    out
  }
  
  # Build the table fields from your simulation output columns
  tab <- data.frame(
    n         = df$ns,
    Pi0       = dgp_label(as.character(df$dgps)),
    UncPwr    = 1 - df$beta_1,
    Delta_c   = df$delta0,
    MeanHat   = df$Mean_deltahat,
    SDHat     = df$SD_deltahat,
    MeanSE    = df$SD_EST_deltahat,
    Cover95   = df$Cover_deltahat,
    stringsAsFactors = FALSE
  )
  
  # Enforce paper ordering: n ascending, then DGP in specified order
  tab$dgps_internal <- as.character(df$dgps)
  tab$dgps_internal <- factor(tab$dgps_internal, levels = dgp_order)
  tab <- tab[order(tab$n, tab$dgps_internal), ]
  tab$dgps_internal <- NULL
  
  # Two-decimal formatting exactly like the paper
  fmt2 <- function(x) sprintf("%.2f", x)
  
  # Compose LaTeX lines
  lines <- c(
    "\\begin{tabular}{ r l r r r r r r r }",
    "\\hline \\hline ",
    "\\multicolumn{4}{c}{Parameters }& & \\multicolumn{4}{c}{Results}\\\\",
    "\\cline{1-4}  \\cline{6-9}\\\\",
    "  $n$ & $\\Pi_0$ & Unc. Pwr. & $\\Delta_c$ & & Mean $\\hat{\\Delta}_{c,n}$ & SD $\\hat{\\Delta}_{c,n}$ & Mean St. Err. & Cover of 95\\% CI\\\\",
    "\\hline "
  )
  
  # Row writer with an empty 5th column to match the paper (`& &`)
  row_tex <- function(r) {
    paste0(
      "  ", r$n, "   &   ", r$Pi0,
      " & ", fmt2(r$UncPwr),
      " & ", fmt2(r$Delta_c),
      " & & ", fmt2(r$MeanHat),
      " & ", fmt2(r$SDHat),
      " & ", fmt2(r$MeanSE),
      " & ", fmt2(r$Cover95),
      "\\\\"
    )
  }
  
  # Write n=50 block then divider, then n=500 block (as in paper)
  ns_sorted <- sort(unique(tab$n))
  if (length(ns_sorted) == 0) stop("No rows to write.")
  if (length(ns_sorted) > 2) {
    # still works, but paper example shows two panels; keep behavior explicit
    warning("More than two n values detected; table will include all panels separated by \\hline.")
  }
  
  for (k in seq_along(ns_sorted)) {
    nn <- ns_sorted[k]
    block <- tab[tab$n == nn, ]
    lines <- c(lines, vapply(seq_len(nrow(block)), function(i) row_tex(block[i, ]), character(1)))
    
    # Panel separator exactly like paper (a single \hline between panels)
    if (k < length(ns_sorted)) {
      lines <- c(lines, "\\hline  ")
    }
  }
  
  lines <- c(lines, "\\hline", "\\end{tabular}")
  
  print(tex_file)
  writeLines(lines, tex_file)
  invisible(tab)
  
}
# --- Application MM -> paper-exact LaTeX table (base R, no kableExtra) ---

write_mm_table_paper_exact <- function(
    out_RCT_tscores, out_NE_tscores, out_RCT_articles, out_NE_articles,
    J_RCT_tscores, J_NE_tscores, J_RCT_articles, J_NE_articles,
    eps_RCT_tscores, eps_NE_tscores, eps_RCT_articles, eps_NE_articles,
    t_RCT, t_NE,
    tex_file,
    cv = 1.96,
    negate_delta = TRUE,
    # Optional zcurve block (paper shows only under "By Number of t-scores")
    zcurve_delta_RCT = NULL, zcurve_se_RCT = NULL,
    zcurve_delta_NE  = NULL, zcurve_se_NE  = NULL,
    digits_delta = 3, digits_se = 3, digits_ci = 2,
    digits_theta = 2, digits_theta_se = 2, digits_power = 2,
    digits_pval = 2
) {
  
  # ---------- helpers ----------
  # format numeric with "paper style" (drop leading 0 for |x|<1 => .072, -.072, .00)
  fmt <- function(x, d) {
    s <- sprintf(paste0("%.", d, "f"), x)
    s <- sub("^(-?)0\\.", "\\1.", s)  # 0.12 -> .12 ; -0.12 -> -.12
    s
  }
  fmt_paren <- function(x, d) paste0("(", fmt(x, d), ")")
  fmt_ci <- function(lo, hi, d) paste0("[", fmt(lo, d), ", ", fmt(hi, d), "]")
  fmt_eps <- function(x, d = 2) paste0("$", fmt(x, d), "$")
  
  # ---------- core quantities from estimator outputs (as in Application MM) ----------
  # Application MM prints results from out_* and uses varest_delta, varest_theta, etc. [1](https://gtvault-my.sharepoint.com/personal/sfaridani6_gatech_edu/Documents/Microsoft%20Copilot%20Chat%20Files/Application%20MM.txt)
  outs <- list(
    RCT_ts  = out_RCT_tscores,
    NE_ts   = out_NE_tscores,
    RCT_art = out_RCT_articles,
    NE_art  = out_NE_articles
  )
  
  deltas <- vapply(outs, function(o) o$deltahat, numeric(1))
  ses    <- sqrt(vapply(outs, function(o) o$varest_delta, numeric(1)))
  ci_lo  <- deltas - 1.96 * ses
  ci_hi  <- deltas + 1.96 * ses
  
  thetas    <- vapply(outs, function(o) o$thetahat, numeric(1))
  theta_ses <- sqrt(vapply(outs, function(o) o$varest_theta, numeric(1)))
  
  # Paper table displays positive Δ (your script often negates when printing “results”). [1](https://gtvault-my.sharepoint.com/personal/sfaridani6_gatech_edu/Documents/Microsoft%20Copilot%20Chat%20Files/Application%20MM.txt)
  if (negate_delta) {
    deltas2 <- -deltas
    ci_lo2  <- -ci_hi
    ci_hi2  <- -ci_lo
    deltas <- deltas2
    ci_lo  <- ci_lo2
    ci_hi  <- ci_hi2
  }
  
  # p-values for equality RCT vs NE under each tuning regime (uses your pval_comp()) [1](https://gtvault-my.sharepoint.com/personal/sfaridani6_gatech_edu/Documents/Microsoft%20Copilot%20Chat%20Files/Application%20MM.txt)[2](https://deepwiki.com/haozhu233/kableExtra/4.2-grouping-rows)
  p_ts  <- pval_comp(out_RCT_tscores, out_NE_tscores)
  p_art <- pval_comp(out_RCT_articles, out_NE_articles)
  
  # Status quo power: empirical rejection rate at cv in observed t-scores (RCT and NE) [1](https://gtvault-my.sharepoint.com/personal/sfaridani6_gatech_edu/Documents/Microsoft%20Copilot%20Chat%20Files/Application%20MM.txt)
  pwr_RCT <- mean(abs(t_RCT) > cv, na.rm = TRUE)
  pwr_NE  <- mean(abs(t_NE)  > cv, na.rm = TRUE)
  
  # counts from estimator outputs (these are computed inside estimator()) [2](https://deepwiki.com/haozhu233/kableExtra/4.2-grouping-rows)[1](https://gtvault-my.sharepoint.com/personal/sfaridani6_gatech_edu/Documents/Microsoft%20Copilot%20Chat%20Files/Application%20MM.txt)
  n_t_RCT <- out_RCT_tscores$num_tscores
  n_t_NE  <- out_NE_tscores$num_tscores
  n_a_RCT <- out_RCT_tscores$num_articles
  n_a_NE  <- out_NE_tscores$num_articles
  
  # ---------- build LaTeX lines EXACTLY in your provided format ----------
  # Note: we write literal LaTeX; your pasted input had HTML &amp; but TeX needs &.
  
  lines <- c(
    "\\begin{tabular}{l c c c c c  }",
    "\\hline \\hline ",
    " & \\multicolumn{2}{c}{\\textit{By Number of $t$-scores}} &  &\\multicolumn{2}{c}{\\textit{By Number of Articles}} \\\\",
    "\\cline{2-3} \\cline{5-6} ",
    "\\\\ & RCTs & Natural Experiments &  & RCTs & Natural Experiments \\\\",
    "\\hline"
  )
  
  # Row block for Delta_hat
  lines <- c(lines,
             paste0(" $\\hat{\\Delta}$ & ", fmt(deltas["RCT_ts"], digits_delta), " & ", fmt(deltas["NE_ts"], digits_delta),
                    "  & & ", fmt(deltas["RCT_art"], digits_delta), " & ", fmt(deltas["NE_art"], digits_delta), "\\\\"),
             paste0("  & ", fmt_paren(ses["RCT_ts"], digits_se), "  & ", fmt_paren(ses["NE_ts"], digits_se),
                    " & & ", fmt_paren(ses["RCT_art"], digits_se), " & ", fmt_paren(ses["NE_art"], digits_se), " \\\\"),
             paste0("  & ", fmt_ci(ci_lo["RCT_ts"], ci_hi["RCT_ts"], digits_ci), " & ", fmt_ci(ci_lo["NE_ts"], ci_hi["NE_ts"], digits_ci),
                    " & & ", fmt_ci(ci_lo["RCT_art"], ci_hi["RCT_art"], digits_ci), " & ", fmt_ci(ci_lo["NE_art"], ci_hi["NE_art"], digits_ci),
                    "\\\\ ")
  )
  
  # Row block for theta_hat
  lines <- c(lines,
             paste0("  $\\hat{\\theta}$ & ", fmt(thetas["RCT_ts"], digits_theta), " & ", fmt(thetas["NE_ts"], digits_theta),
                    " & & ", fmt(thetas["RCT_art"], digits_theta), " & ", fmt(thetas["NE_art"], digits_theta), "\\\\"),
             paste0("    & ", fmt_paren(theta_ses["RCT_ts"], digits_theta_se), " & ", fmt_paren(theta_ses["NE_ts"], digits_theta_se),
                    " & & ", fmt_paren(theta_ses["RCT_art"], digits_theta_se), " & ", fmt_paren(theta_ses["NE_art"], digits_theta_se), " \\\\ ")
  )
  
  lines <- c(lines,
             " \\hline  ",
             paste0("  J  & ", as.integer(round(J_RCT_tscores)), "& ", as.integer(round(J_NE_tscores)),
                    " &  & ", as.integer(round(J_RCT_articles)), " & ", as.integer(round(J_NE_articles)), "\\\\"),
             paste0("  $\\epsilon$ & ", fmt_eps(eps_RCT_tscores, 2), " & ", fmt_eps(eps_NE_tscores, 2),
                    " &  & ", fmt_eps(eps_RCT_articles, 2), " & ", fmt_eps(eps_NE_articles, 2), " \\\\ \\hline"),
             paste0("     $=$ RCT (p-value) &  & ", fmt(p_ts, digits_pval),
                    " &  &  &", fmt(p_art, digits_pval), " \\\\ \\hline ")
  )
  
  # zcurve block (paper shows values only under t-score tuning; article columns blank)
  if (!is.null(zcurve_delta_RCT) && !is.null(zcurve_se_RCT) &&
      !is.null(zcurve_delta_NE)  && !is.null(zcurve_se_NE)) {
    
    lines <- c(lines,
               paste0("      $\\hat{\\Delta}$ \\verb|zcurve| & ", fmt(zcurve_delta_RCT, digits_delta),
                      " &  ", fmt(zcurve_delta_NE, digits_delta), " & &  & \\\\"),
               paste0("  & ", fmt_paren(zcurve_se_RCT, digits_se),
                      "  & ", fmt_paren(zcurve_se_NE, digits_se), " & &  &  \\\\ \\hline ")
    )
    
  } else {
    # exact blank structure as your template
    lines <- c(lines,
               "      $\\hat{\\Delta}$ \\verb|zcurve| &  &   & &  & \\\\",
               "  &   &  & &  &  \\\\ \\hline "
    )
  }
  
  # bottom block
  lines <- c(lines,
             paste0("        Status Quo Power & ", fmt(pwr_RCT, digits_power), " & ", fmt(pwr_NE, digits_power),
                    " && ", fmt(pwr_RCT, digits_power), " & ", fmt(pwr_NE, digits_power), " \\\\"),
             paste0("  $t$-scores & ", as.integer(round(n_t_RCT)), " & ", as.integer(round(n_t_NE)),
                    " &  & ", as.integer(round(n_t_RCT)), " & ", as.integer(round(n_t_NE)), "  \\\\"),
             paste0("  Articles & ", as.integer(round(n_a_RCT)), "& ", as.integer(round(n_a_NE)),
                    " &  & ", as.integer(round(n_a_RCT)), "& ", as.integer(round(n_a_NE)), " \\\\"),
             "\\hline",
             "\\end{tabular}"
  )
  
  # ---------- robust save ----------
  dir.create(dirname(tex_file), recursive = TRUE, showWarnings = FALSE)
  con <- file(tex_file, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con = con)
  message("MM paper table written to: ", normalizePath(tex_file, winslash = "/", mustWork = FALSE))
  
  invisible(lines)
}


write_ml_table_paper_exact <- function(
    # estimator() outputs
  ML_by_tscores, ML_by_tscores_r,
  ML_by_sites,   ML_by_sites_r,
  
  # tuning params (already computed in your script)
  J_ML_tscores,   J_ML_tscores_r,
  J_ML_sites,     J_ML_sites_r,
  eps_ML_tscores, eps_ML_tscores_r,
  eps_ML_sites,   eps_ML_sites_r,
  
  # constants shown in the table
  D_main  , D_rob  ,
  C_main ,    C_rob ,
  
  # data used for unconditional power + counts
  tscores, sites, treatments,
  
  # output
  tex_file,
  cv = 1.96,
  
  # formatting
  digits_delta = 3,
  digits_se    = 3,
  digits_ci    = 3,
  digits_theta = 2,
  digits_theta_se = 2,
  digits_power = 2,
  
  # paper style: CI lower bound truncated at 0
  truncate_ci_at_zero = TRUE
) {
  # ---- helpers: paper-style numeric formatting ----
  fmt <- function(x, d) {
    s <- sprintf(paste0("%.", d, "f"), x)
    s <- sub("^(-?)0\\.", "\\1.", s)  # 0.081 -> .081 ; -0.081 -> -.081
    s
  }
  fmt_paren <- function(x, d) paste0("(", fmt(x, d), ")")
  
  fmt_ci <- function(lo, hi, d) {
    # match paper: show "0" (not .000) when truncated to zero
    lo_str <- ifelse(abs(lo) < 1e-12, "0", fmt(lo, d))
    paste0("[", lo_str, ", ", fmt(hi, d), "]")
  }
  
  fmt_eps <- function(x, d = 2) paste0("$", fmt(x, d), "$")
  
  # ---- extract results from estimator outputs (as in Application ML) ----
  outs <- list(
    ts_main = ML_by_tscores,
    ts_rob  = ML_by_tscores_r,
    st_main = ML_by_sites,
    st_rob  = ML_by_sites_r
  )
  
  deltas <- vapply(outs, function(o) o$deltahat, numeric(1))
  ses    <- sqrt(vapply(outs, function(o) o$varest_delta, numeric(1)))
  ci_lo  <- deltas - 1.96 * ses
  ci_hi  <- deltas + 1.96 * ses
  
  if (truncate_ci_at_zero) {
    ci_lo <- pmax( ci_lo,0)
  }
  
  thetas    <- vapply(outs, function(o) o$thetahat, numeric(1))
  theta_ses <- sqrt(vapply(outs, function(o) o$varest_theta, numeric(1)))
  
  # ---- bottom-block quantities ----
  unc_power <- mean(abs(tscores) > cv, na.rm = TRUE)
  
  n_tscores   <- length(tscores)
  n_sites     <- length(unique(sites))
  n_treat     <- length(unique(treatments))
  
  # ---- table parameters ----
  Js <- c(ts_main = J_ML_tscores, ts_rob = J_ML_tscores_r,
          st_main = J_ML_sites,   st_rob = J_ML_sites_r)
  
  eps <- c(ts_main = eps_ML_tscores, ts_rob = eps_ML_tscores_r,
           st_main = eps_ML_sites,   st_rob = eps_ML_sites_r)
  
  # ---- IMPORTANT: your template declares 7 columns: {l c c c c c c}
  # To avoid LaTeX alignment errors, we write an EXTRA blank final column on every row.
  # That means every row ends with "... & \\\\" (i.e., a trailing empty cell).
  endrow <- " \\\\"
  
  # ---- build LaTeX lines EXACTLY like the provided format ----
  lines <- c(
    "\\begin{tabular}{l c c c c c c }",
    "\\hline \\hline ",
    " & \\multicolumn{2}{c}{\\textit{By Number of $t$-scores}} &  &\\multicolumn{2}{c}{\\textit{By Number of Sites}} & \\\\",
    "\\cline{2-3} \\cline{5-6} ",
    "\\\\ & Main Specification & Rob. Check &  & Main Specification & Rob. Check & \\\\",
    "\\hline"
  )
  
  # Δ̂c block
  lines <- c(lines,
             paste0(" $\\hat{\\Delta}_c$ & ",
                    fmt(deltas["ts_main"], digits_delta), " & ", fmt(deltas["ts_rob"], digits_delta),
                    " & & ",
                    fmt(deltas["st_main"], digits_delta), " & ", fmt(deltas["st_rob"], digits_delta),
                    " &", endrow),
             
             paste0("& ", fmt_paren(ses["ts_main"], digits_se), " & ", fmt_paren(ses["ts_rob"], digits_se),
                    " & & ",
                    fmt_paren(ses["st_main"], digits_se), " & ", fmt_paren(ses["st_rob"], digits_se),
                    " &", endrow),
             
             paste0(" &  ", fmt_ci(ci_lo["ts_main"], ci_hi["ts_main"], digits_ci),
                    " & ", fmt_ci(ci_lo["ts_rob"],  ci_hi["ts_rob"],  digits_ci),
                    " & & ",
                    fmt_ci(ci_lo["st_main"], ci_hi["st_main"], digits_ci),
                    " & ", fmt_ci(ci_lo["st_rob"],  ci_hi["st_rob"],  digits_ci),
                    " &", endrow),
             
             "    \\hline "
  )
  
  # θ̂ block
  lines <- c(lines,
             paste0("  $\\hat{\\theta}$ & ",
                    fmt(thetas["ts_main"], digits_theta), " & ", fmt(thetas["ts_rob"], digits_theta),
                    " & & ",
                    fmt(thetas["st_main"], digits_theta), " & ", fmt(thetas["st_rob"], digits_theta),
                    " &", endrow),
             
             paste0(" & ", fmt_paren(theta_ses["ts_main"], digits_theta_se), " & ", fmt_paren(theta_ses["ts_rob"], digits_theta_se),
                    " & & ",
                    fmt_paren(theta_ses["st_main"], digits_theta_se), " & ", fmt_paren(theta_ses["st_rob"], digits_theta_se),
                    " &", endrow),
             
             " \\hline  "
  )
  
  # D, C, J, eps rows (exact labels/structure)
  lines <- c(lines,
             paste0(" D  &  ", format(D_main, scientific = TRUE), " &  ", format(D_rob, scientific = TRUE),
                    " & & ", format(D_main, scientific = TRUE), " & ", format(D_rob, scientific = TRUE),
                    " &", endrow),
             
             paste0("  C  &  $", C_main, "$ & $", C_rob, "$ & & $", C_main, "$ & $", C_rob, "$",
                    " &", endrow),
             
             paste0("  J    & ", as.integer(round(Js["ts_main"])), " & ", as.integer(round(Js["ts_rob"])),
                    " & & ", as.integer(round(Js["st_main"])), "& ", as.integer(round(Js["st_rob"])),
                    " &", endrow),
             
             paste0("  $\\epsilon$  & ", fmt_eps(eps["ts_main"], 2), " & ", fmt_eps(eps["ts_rob"], 2),
                    " & &  ", fmt_eps(eps["st_main"], 2), " & ", fmt_eps(eps["st_rob"], 2),
                    " &", endrow),
             
             " \\hline"
  )
  
  # Bottom block: unconditional power + counts (repeated across columns)
  lines <- c(lines,
             paste0("     Unconditional Power & ", fmt(unc_power, digits_power), " & ", fmt(unc_power, digits_power),
                    " &  & ", fmt(unc_power, digits_power), " &", fmt(unc_power, digits_power),
                    " &", endrow),
             
             paste0("  $t$-scores & ", n_tscores, " & ", n_tscores, " & & ", n_tscores, " & ", n_tscores,
                    " &", endrow),
             
             paste0("  Sites & ", n_sites, "& ", n_sites, " & & ", n_sites, " & ", n_sites,
                    " &", endrow),
             
             paste0("  Treatments & ", n_treat, " & ", n_treat, "& & ", n_treat, " & ", n_treat,
                    " &", endrow),
             
             "\\hline",
             "\\end{tabular}"
  )
  
  # ---- robust write ----
  dir.create(dirname(tex_file), recursive = TRUE, showWarnings = FALSE)
  con <- file(tex_file, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con = con)
  message("ML paper table written to: ", normalizePath(tex_file, winslash = "/", mustWork = FALSE))
  
  invisible(lines)
}