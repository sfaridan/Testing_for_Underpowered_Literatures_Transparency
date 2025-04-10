library(Rsolnp)

ngrid <- 25
tol <- .001
tgrid <-grid(8, ngrid)
hgrid <- grid(5,ngrid)
tgrid <- tgrid[tgrid >= 0]
hgrid <- hgrid[hgrid >= 0]
V <- cdf_mat(tgrid,hgrid,c=1,sided=2)
Vc <- cdf_mat(tgrid,hgrid,c=(2),sided=2)
svdV <- svd(V)
eigenleft <- sum((abs(svdV$d))*(abs(svdV$d) < tol) )
maxbias <- 1/(sqrt(2*pi)*ngrid)+eigenleft
maxbias
values_inv <- rep(0,length(svdV$d))
values_inv[abs(svdV$d)>tol] <- 1/Re(svdV$d[abs(svdV$d)>tol])
values_inv <- 1/(1+Re(svdV$d)) 
Lambda_inv <- diag(values_inv)
vals_reg <- svdV$d * (abs(svdV$d) > tol) 
Vreg <- svdV$u%*%( diag(vals_reg) )%*%t(svdV$v)

pis<-rep(0,length(hgrid))
pis[sum(hgrid<=1.5)]<-1/3
pis[sum(hgrid<=1.9)]<-1/3
pis[sum(hgrid<=2.1)]<-1/3
Vp <- V%*%pis
Vcp <- Vc%*%pis
invmat <- svdV$v %*% Lambda_inv %*% t(svdV$u)
pihat <- invmat %*% V %*%pis
plot(tgrid,Vp,type="l")
lines(tgrid,V%*%pihat,col="red")
plot(tgrid,Vcp,type="l")
lines(tgrid,Vc%*%pihat,col="red")
c(norm(Vc%*%(pis-pihat)),maxbias)

Tscore_data <- c(rnorm(100,0,1),rnorm(100,2,1)) 
ecdf <- empcdf_bygrid(Tscore_data,tgrid)
pihat_noise <- invmat %*% ecdf
plot(tgrid,Vp,type="l")
lines(tgrid,ecdf,col="red")
plot(tgrid,Vp,type="l")
lines(tgrid,V%*%pihat_noise,col="red")
plot(tgrid,Vcp,type="l")
lines(tgrid,Vc%*%pihat_noise,col="red")

#asymptotic normality?
nsims <- 1000
ests <- rep(0,nsims)
ests2 <- rep(0,nsims)
zscores <- rep(0,nsims)
tc <-  1.96
index <- sum(tgrid<=tc)
estmat <- Vc[index,]%*%invmat
truetheta <- Vc[index,]%*%pis
for (i in 1:nsims){
  Tscore_data <- c(rnorm(100,1.9,1),rnorm(100,1.5,1),rnorm(100,2.1,1)) 
  rand <- runif(length(Tscore_data))
  #Tscore_data <- Tscore_data[abs(Tscore_data) >= rep(tc,length(Tscore_data)) | rand <= 0.75 ]
  Tscore_data <- abs(Tscore_data)
  ecdf <- empcdf_bygrid(Tscore_data,tgrid)
  obj<- function(pihat){
    return(norm(V%*%pihat-ecdf))
  }
  objpb<- function(pihat){
    
    ecdf_pb <- pub_bias(tgrid, ecdf,tc,1/pihat[length(pihat)])
    penalty <- (ecdf_pb[sum(tgrid <= tc+0.05)] +ecdf_pb[sum(tgrid <= tc-0.05)] - 2*ecdf_pb[sum(tgrid <= tc)])/.05 
    return(norm(V%*%pihat[1:length(hgrid)]-ecdf_pb)+0.001*abs(penalty)  )
  }
  objtrue<- function(pihat){
    tcdf <-  pub_bias(tgrid, V%*%pis,tc,0.75)
    tcdf_pb <- pub_bias(tgrid, tcdf,tc,1/pihat[length(pihat)])
    penalty <- (tcdf_pb[sum(tgrid <= tc+0.05)] +tcdf_pb[sum(tgrid <= tc-0.05)] - 2*tcdf_pb[sum(tgrid <= tc)])/.05 
    return(norm(V%*%pihat[1:length(hgrid)]-tcdf_pb)+0.001*abs(penalty)  )
  }
  ests[i] <- estmat%*%ecdf
  X0 <- c(rep(1/length(hgrid),length(hgrid)),0.9)
  X0 <- c(rep(1/length(hgrid),length(hgrid)))
  sols <- solnp(X0, obj,LB = c(rep(0,length(hgrid))),control = list(trace = 0))
  ests2[i] <-Vc[index,]%*%sols$pars[1:length(hgrid)]
  #ests2[i] <- Vc[index,]%*%invmat%*%V%*%sols$pars[1:length(hgrid)]
  zscores[1:i] <- (ests[1:i]-mean(ests[1:i]))/sd(ests[1:i])
  print(c(i,nsims))
  print(c(sols$pars[length(hgrid)+1],mean(ests[1:i]),truetheta, sd(ests[1:i]),mean(abs(zscores[1:i])<= 1.65),mean(abs(zscores[1:i])<= 1.96),mean(abs(zscores[1:i])<= 2.58) ))
  #print(c(sols$pars[length(hgrid)+1],mean(ests2[1:i]),truetheta, sd(ests2[1:i]),mean(abs(zscores[1:i])<= 1.65),mean(abs(zscores[1:i])<= 1.96),mean(abs(zscores[1:i])<= 2.58)
  print(c(mean(ests[1:i]),mean(ests2[1:i]),truetheta,sd(ests[1:i]),sd(ests2[1:i]) ))
}
hist(zscores)
#0.01044277/sqrt(10)
#[1] 0.003302294
print(c(sd(ests)*1.96,mean(ests)-truetheta,maxbias))

estimator_discretized <- function(Tscore_data, t, cs,sided=1,tol=10e-6,tupper= inf){
  tgrid  <- grid(max(abs(Tscore_data))+4,10)
  V <- cdf_mat(tgrid,tgrid,c=1,sided=1)
  Vc <- cdf_mat_overcs(t,tgrid,cs,sided=2) 
  
  eigenV <- eigen(V)
  
  vals <- eigenV$values * (abs(eigenV$values) >=  tol)
  Lambda <- diag(vals)
  
}

estimator <- function(Tscore_data,articles,tc,cs,lambda=1.452024,pb=1,penalty=0,sided=1,nsims_inf=100,tupper = inf){
  t <- tc
  Maxt <- max(Tscore_data) 
  Mt    <- min(c(tupper,Maxt))
  hgrid  <- c(grid(Mt+4,10),(round(Mt+4):round(Mt+50)  ))
  tgrid <- c(grid(Mt,120))
  if (sided == 2){
    Tscore_data <- abs(Tscore_data)
    hgrid <- hgrid[ hgrid >= 0 ]
    tgrid <- tgrid[ tgrid >= 0 ]
  }
  V <- cdf_mat(tgrid,hgrid,c=1,sided=2)
  Vc <- cdf_mat_overcs(t,hgrid,cs,sided=2) 
  V1 <- cdf_mat_overcs(t,hgrid,rep(1,length(cs)),sided=2)
  out <- estfit_cdf(Tscore_data,tgrid,V,pb=pb,tc=tc,sided=sided,penalty=penalty)
  pihats <- out$pi_hat
  out$ests <- Vc%*%pihats
  out$U <- Utility(out$ests,cs,lambda=lambda)
  out$V <- V
  ci_sims<- resample_est(Tscore_data,articles,tgrid,V,out$pi_hat,nsims_inf,pb=pb,tc=tc,sided=sided,penalty=penalty)
  cis <- est_ci(ci_sims, Vc,V1,pb,out$raw, cs,lambda=lambda, size=0.05)
  
  out$alpha_hat <- out$raw[length(out$raw)]
  out$cs <- cs
  out$ugain <- cis$ugain
  out$citop <- cis$citop
  out$cibot <- cis$cibot
  out$utop <- cis$utop
  out$ubot <- cis$ubot
  out$alpha_top <- cis$alpha_top
  out$alpha_bot <- cis$alpha_bot
  out$tgrid <- tgrid
  out$hgrid <- hgrid
  out$pb <- pb
  out$tc <- tc
  out$t <- t
  return(out)
}

Utility <- function(ests,cs,lambda=1.452024){
  return( ((1+lambda)*(1-ests)-1 )/cs  )
}

#construct a matrix with columns corresponding to t-scores in t-grid
# columns corresponding to h values in hgrid and c is the multiplier on h
# output: V%*%pi is the CDF of the t-score
cdf_mat <-  function(tgrid,hgrid,c=1,sided=1){
  
  argst = matrix(rep(tgrid,length(hgrid)),ncol=length(hgrid))
  argsh = t(matrix(rep(hgrid,length(tgrid)),ncol=length(tgrid)))

  V = pnorm(argst - sqrt(c)*argsh)

  if (sided == 2){
    V = V- pnorm(-argst - sqrt(c)*argsh)
  }
  
  return(V)
}

cdf_mat_overcs <-  function(t,hgrid,cs,sided = 1){
  
  argscs = matrix(rep(cs,length(hgrid)),ncol=length(hgrid))
  argsh = t(matrix(rep(hgrid,length(cs)),ncol=length(cs)))
  
  V = pnorm(t - sqrt(argscs)*argsh)
  if (sided == 2){
    V = V- pnorm(-t - sqrt(argscs)*argsh)
  }
  
  return(V)
}

#Finds pi closest to the empirical CDF of the data
# if pub bias is included, the last element of the solution vector is the pb parameter
estfit_cdf <- function(Tscore_data,tgrid,V,pb=0,tc=1.65,sided=1,penalty=0){
  n <- length(Tscore_data)
  hs <- dim(V)[2]
  if (sided==2){
    Tscore_data <-  abs(Tscore_data) 
  }
  
  #compute empirical cdf
  ecdf_data<- empcdf_bygrid(Tscore_data,tgrid)

  
  #Set up the minimization problem
  X0 <- rep(1/hs,hs) #initial guess of distribution pi
  LB <- rep(0,hs) #lower bound for pi
  UB <- rep(1,hs) #upper bound for pi
  obj <- function(x){
    return( norm( V%*%x - ecdf_data ,"2" )  )
  }
  eqfun <- function(x){
    return(sum(x[1:hs]))
  }
  eqB <- 1
  
  #If we are including pub bias
  if (pb==1){
    index <- sum(tgrid<=tc)
    X0 <- c(X0,0.5) #add pb parameter
    LB <- c(LB,0.01)
    UB <- c(UB,1) #ub of 1 makes sense but instead use 2 to alert us to weird solutions
    
    obj <- function(x){
      tcdf_hat <- pub_bias(tgrid, ecdf_data,tc,1/x[hs+1])
      
      fit <- norm( V%*%(x[1:hs]) - tcdf_hat ,"2" )
      pen <- 1*penalty*abs( (tcdf_hat[index]-tcdf_hat[index-1]) -(tcdf_hat[index+1]-tcdf_hat[index]) )

      return( fit +pen )
    }
  }
  
  #solve problem

  xsol <- solnp(X0, obj,eqfun=eqfun,eqB=eqB, LB = LB, UB =  UB,control = list(trace = 0))

  
  #prepare output
  output <- list()
  output$raw <- xsol$pars
  output$pi_hat <- xsol$pars[1:hs]
    
  return(output)
}


resample_est  <- function(Tscore_data,articles,tgrid,V,ests,nsims_inf=100,pb=0,tc=1.65,sided=1,penalty=0){
  dimsol <- dim(V)[2]+pb #elements in solution
  sol_sims <- matrix(rep(0,dimsol*nsims_inf),nrow = dimsol)
  unique_articles <- unique(articles)
  framd <- data.frame(Tscore_data,articles)
   for (sim in 1:nsims_inf){
     #
    framd_resampled <- do.call(rbind, sample(split(framd, framd$articles), replace = T))
    resampled_data <- framd_resampled$Tscore_data
    out <- estfit_cdf(resampled_data,tgrid,V,pb,tc,penalty=penalty)
    sol_sims[,sim] <- out$raw-ests

  }
  return(sol_sims)
}

# given simulations returns ci covering counterfactual cdf 
# with factor c
est_ci <- function(ci_sims, Vc,V1,pb, outs, cs, lambda=1.452024, size=0.05){
  Vcx <- Vc-V1
  cis <- matrix(rep(0,2*(dim(ci_sims)[1])),ncol=2)
  # columns are values of h, rows are values of c
  
  nsims <- dim(ci_sims)[2]
  ncs <- dim(Vc)[2]
  pts <- Vcx%*%ci_sims[1:ncs,] 
  alphas <- ci_sims[dim(ci_sims)[1],] 
  ptests <- Vc%*%(ci_sims[1:ncs,]+outs[1:ncs])
  cmaximizer <- rep(0,nsims)
  pmaximizer <- rep(0,nsims)
  c1 <- sum(cs<=1)
  Ugainsim <- matrix(rep(0,nsims*length(cs)),nrow=length(cs))
  for (sim in 1:nsims)  
  {
    U <- Utility(ptests[,sim],cs,lambda=lambda)
    cmaximizer[sim]  <-  cs[which.max(U)]
    pmaximizer[sim]  <-  1-ptests[which.max(U),sim]
    Ugainsim[,sim] <- U -U[c1]
  }

  
  qtop <- function(x){
    return(quantile(x, 1-size/2))
  }
  qbot <- function(x){
    return(quantile(x, size/2))
  }
  
  tops <- abs(apply(pts, 1, qtop))
  bots <- abs(apply(pts, 1, qbot))
  utops <- abs(qtop(cmaximizer))
  ubots <- abs(qbot(cmaximizer))
  ptops <- abs(qtop(pmaximizer))
  pbots <- abs(qbot(pmaximizer))
  atops <- abs(qtop(alphas))
  abots <- abs(qbot(alphas))

  
  

  
  ests <- Vcx%*%outs[1:dim(Vc)[2]]
  Uestimated <- Utility(Vc%*%outs[1:ncs],cs,lambda=lambda)
  cmaxest <- cs[which.max(Uestimated)]
  pmaxest <- 1-(Vc%*%outs[1:ncs])[which.max(Uestimated)]
  alpha <- outs[length(outs)]
  Ugain <- Uestimated - Utility((V1%*%outs[1:ncs])[1],1,lambda=lambda)
  

  
  #symmeterize the cis
  delta <- apply( rbind(tops,bots) ,2,max  )
  udelta <- max(c(abs(utops -cmaxest),abs(ubots -cmaxest)))
  pdelta <- max(c(abs(ptops -pmaxest),abs(pbots -pmaxest)))
  adelta <- max(c(abs(atops -alpha),abs(abots -alpha)))
  
  udeltatop <- abs(utops-cmaxest)
  udeltabot <- abs(ubots-cmaxest)
  
  cis$citop <- ests + delta
  cis$cibot <- ests - delta
  cis$utop <- cmaxest+max(udeltabot,udeltatop) #utops #cmaxest + udelta
  cis$ubot <- cmaxest - max(udeltabot,udeltatop)  #ubots #cmaxest - udelta
  cis$ptop <- pmaxest+pdelta
  cis$pbot <- pmaxest-pdelta
  cis$alpha_top <- alpha + adelta
  cis$alpha_bot <- alpha - adelta
  
  cis$ugain <-  apply(1*(Ugainsim>0),1,mean )

  
  
  return(cis)
}

empcdf_bygrid<- function(Tscore_data,tgrid){
  
  ecdf <- 0*tgrid
  for (t in 1:length(ecdf)){
    ecdf[t] = mean(Tscore_data<= tgrid[t])
  }
  return(ecdf)
}

pub_bias <- function(tgrid, cdf_prepb,tc,pb_param){
  ts <- length(tgrid)
  index <- sum(tgrid<=tc)
  pre <- cdf_prepb[index]
  pmf_prepb <- c(cdf_prepb[1],cdf_prepb[2:ts]-cdf_prepb[1:(ts-1)] )
  pmf_with <- pmf_prepb*((tgrid<=tc)*pb_param+(tgrid>tc)) / (pb_param*pre+(1-pre)   )
  cdf_withpb <- cumsum(pmf_with)
}

# compute pihat and cis for several cs
sim_param <- function(tgrid,hgrid,pis,pb_param,tc,n,nsims,pb,t,cs,lambda=1.452024,sided=1,penalty=0,nsims_inf=100){
  V <- cdf_mat(tgrid,hgrid,1,sided=sided)
  cdf_true <-  V%*%pis
  cdf_pb <- pub_bias(tgrid, cdf_true,tc, pb_param)
  plot(tgrid, cdf_pb,type="l")
  plot(tgrid,pub_bias(tgrid,cdf_pb,tc,1/pb_param),type="l")
  index <- sum(tgrid<=t)
  last <- length(hgrid)+1
  pihats <- matrix(rep(0,nsims*( length(hgrid)+ pb )),nrow=nsims)
  covers_alpha <- rep(0,nsims)
  articles <- 1:n #expand on this later
  Vc <- cdf_mat_overcs(t,hgrid,(cs),sided=sided) 
  V1 <- cdf_mat_overcs(t,hgrid,rep(1,length(cs)),sided=sided) 
  truevals <- (Vc-V1)%*%pis
  Utils <- Utility(Vc%*%pis,cs,lambda=lambda)
 
  optim <- cs[which.max(Utils)]
  poptim <- 1-(Vc%*%pis)[which.max(Utils)]

  optim_hat <- rep(0,nsims)
  poptim_hat <- rep(0,nsims)
  covers <- matrix(rep(0,nsims*dim(Vc)[1]),nrow=nsims)
  covers_optim <- rep(0,nsims)
  covers_poptim <- rep(0,nsims)
  rej_inc_bad <- matrix(rep(0,nsims*length(cs)),nrow=length(cs) )
  c1 <- sum(cs<=1)
    for (sim in 1:nsims){
      
    if (20*round(sim/20)==sim){  
    print(paste0("Sim ", sim , " of ", nsims))
    }
    Tscore_data <- gensamp(cdf_pb,tgrid,n)
    
    out <- estfit_cdf(Tscore_data,tgrid,V,pb=pb,tc=tc,sided=sided,penalty=penalty)
    pihats[sim,] <- out$raw 
    ci_sims<- resample_est(Tscore_data,articles,tgrid,V, out$pi_hat,nsims_inf=nsims_inf,pb=pb,tc=tc,sided=sided,penalty=penalty)
    cis <- est_ci(ci_sims, Vc,V1, pb, out$raw, cs, lambda=lambda, size=0.05)
      covers[sim,]<- (cis$cibot<=truevals ) & (cis$citop>=truevals)
    covers_optim[sim]<- (cis$ubot<=optim ) & (cis$utop>=optim)
    covers_poptim[sim]<- (cis$pbot<=poptim ) & (cis$ptop>=poptim)
    
    
    Utils_hat <- Utility(Vc%*%out$pi_hat,cs,lambda=lambda)
    optim_hat[sim] <- cs[which.max(Utils_hat)]
    poptim_hat[sim] <- 1-(Vc%*%out$pi_hat)[which.max(Utils_hat)]

    #print(paste0("Pwr: ", round(1-(V1%*%pis)[1],2), " Pwr hat: ",round(1-(V1%*%out$pi_hat)[1] ,2)))
    #print(paste0("C Optim: ", round(optim,2), " c Est: ", round(optim_hat[sim],2), " CI: ",round(cis$ubot,2), ", ", round(cis$utop,2)) )
    #print(paste0("P Optim: ", round(poptim,2), " pEst: ", round(poptim_hat[sim],2), " CI: ",round(cis$pbot,2), ", ", round(cis$ptop,2)) )
    #print(cbind(cs,1-Vc%*%pis,1-Vc%*%out$pi_hat,Utils-Utils[c1],Utils_hat-Utils_hat[c1], cis$ugain))
    
    rej_inc_bad[,sim] <- 1*(cis$ugain >= 0.95)

  }
  
  coverage <- apply(covers, 2, mean)

  
  outs <- list()
  outs$rmse <- sqrt(apply( (Vc%*%t(pihats[,1:length(hgrid)] -t(matrix(rep(pis,nsims),nrow=length(pis))) ))^2, 1, mean ))
  outs$rmse_diff <- sqrt(apply( ((V1-Vc)%*%t(pihats[,1:length(hgrid)] -t(matrix(rep(pis,nsims),nrow=length(pis))) ))^2, 1, mean ))
  
  hist(Vc[2,]%*%t(pihats[,1:length(hgrid)] -t(matrix(rep(pis,nsims),nrow=length(pis)))))
  print((Vc%*%t(pihats[,1:length(hgrid)] -t(matrix(rep(pis,nsims),nrow=length(pis))))))
    
  outs$cover <- coverage
  outs$cover_optim <- mean(covers_optim)
  outs$cover_poptim <- mean(covers_poptim)
  
  outs$rejects_inc_bad <- apply(rej_inc_bad,1,mean)

  outs$increasesGood <- 1*(Utils -Utils[c1] > 0)

  
  print("cs, Increases good?, Reject incrase bad")
  print(cbind(cs,  outs$increasesGood, outs$rejects_inc_bad))
  
  outs$optim <- optim
  outs$poptim <- poptim
  outs$rmse_optim <- sqrt(mean( (optim_hat-optim)^2 )  )
  outs$rmse_poptim <- sqrt(mean( (poptim_hat-poptim)^2 )  )
  outs$estimands <- Vc%*%pis
  if (pb==1){
    outs$rmse_alpha <- sqrt(  mean( (pihats[,length(hgrid)+1]-pb_param)^2  ) )
  }
  plot(tgrid, V%*%(pihats[1,1:length(hgrid)]))
  lines(tgrid, V%*%(pis))
    return(outs)
}

run_sims <- function(parms,cs,t){
  
  nps <- nrow(parms)
  parms_out <- parms
  parms_out$rmse_alpha <- rep(0,nps)
  for (ind in 1:length(cs)){
    c <- cs[ind]
    name_rmse <- paste0("rmse_t",t,"_c",c)
    parms_out[[name_rmse]] <- rep(0,nps)
    name_rmse_diff <- paste0("rmse_diff_t",t,"_c",c)
    parms_out[[name_rmse_diff]] <- rep(0,nps)
    name_cover <- paste0("cover_",c)
    parms_out[[name_cover]] <- rep(0,nps)
    name_estimand <- paste0("estimand_",c)
    parms_out[[name_estimand]] <- rep(0,nps)
  }
  
  #Loop over parameterizations
  for (parm in 1:nps) {
    print(paste0("Parm ", parm , " of ", nps))
        #Discretization
    tgrid <- grid( parms$Mt[parm], parms$divt[parm])
    hgrid <- grid( parms$Mh[parm], parms$divh[parm])
    if (sided == 2){
      hgrid <- hgrid[ hgrid >= 0 ]
      tgrid <- tgrid[ tgrid >= 0 ]
     
    }
    
    pi_true <- getpi(hgrid,parms$pi_dgps[parm],sided=sided)
    #print(pi_true)
    out <- sim_param(tgrid,hgrid,pi_true,parms$pb_param[parm],parms$tc[parm],parms$N[parm],parms$nsims[parm],parms$pb[parm],t,cs,lambda=parms$lambda[parm],nsims_inf=parms$nsims_inf[parm],sided=parms$sided[parm],penalty=parms$penalty[parm])
    #pihats <- out$pihats
    
    #record results
    for (ind in 1:length(cs)){
      c <- cs[ind]
      name_rmse <- paste0("rmse_t",t,"_c",c)
      rmses_in <- parms_out[[name_rmse]] 
      rmses_in[parm] <-out$rmse[ind]
      parms_out[[name_rmse]] <- rmses_in
      
      name_rmse_diff <- paste0("rmse_diff_t",t,"_c",c)
      rmses_in_diff <- parms_out[[name_rmse_diff]] 
      rmses_in_diff[parm] <-out$rmse_diff[ind]
      parms_out[[name_rmse_diff]] <- rmses_in_diff
      
      name_cover <- paste0("cover_t",t,"_c",c)
      in_cover <- parms_out[[name_cover]] 
      in_cover[parm] <-out$cover[ind]
      parms_out[[name_cover]] <- in_cover
      
      #name_cover <- paste0("cover_",c)
      #parms_out[[name_cover]] <- c(parms_out[[name_cover]],out$cover[ind])
      name_estimand <- paste0("estimand_",c)
      estimand_in <- parms_out[[name_estimand]] 
      estimand_in[parm] <- out$estimands[ind]
      parms_out[[name_estimand]] <- estimand_in
    }
    if (parms$pb==1){
    parms_out$rmse_alpha[parm] <- out$rmse_alpha
    }
    parms_out$cover_optim[parm] <- out$cover_optim
    parms_out$optim[parm] <- out$optim
    parms_out$rmse_optim[parm] <- out$rmse_optim
    parms_out$cover_poptim[parm] <- out$cover_poptim
    parms_out$poptim[parm] <- out$poptim
    parms_out$rmse_poptim[parm] <- out$rmse_poptim
    parms_out$increasesGood[parm] <- out$increasesGood
    parms_out$rejects_inc_bad[parm] <- out$rejects_inc_bad
      
  }
  return(parms_out)
}

grid <- function(M,div){
  g <- ((-M*div):(M*div))/div
  return(g)
}

save_results <- function(parms_out, results_path){
  sanitize_date <- gsub(as.character(Sys.Date()), pattern = "-", replacement = "_")
  sanitize_parms <- paste0("Results", "_from_",length(parms_out$N),"_parms_")
  sanitize_save_name <- paste0(sanitize_parms,sanitize_date)
  path <-paste0(results_path, sanitize_save_name, ".csv") 
  print(path )
  write.csv(x = parms_out, file = path)
}

getpi <- function(hgrid, pi_dgp,sided=1){
  
  index0 = sum(hgrid<=0)
  index1 = sum(hgrid<=1)
  index25 = sum(hgrid<=2.5)
  index28 = sum(hgrid<=2.8)
  index35 = sum(hgrid<=3.5)
  index2 = sum(hgrid<=2)
  pi_true <- 0*hgrid
  
  #half zero, half h=1
  if (pi_dgp==1){
     pi_true[index0] <- 0.5
     pi_true[index25] <- 0.5
  }
  if (pi_dgp==2){
    pi_true[abs(hgrid)<=2.5] <- 1
    pi_true <- pi_true / sum(pi_true)
  }
  if (pi_dgp==3){
    pi_true[index0] <- 1 
  }
  if (pi_dgp==4){
    pi_true[index28] <- 1 
  }
  if (pi_dgp==5){
    pi_true[index2] <- 1 
  }
  if (pi_dgp==6){
    pi_true[index35] <- 1
    #pi_true[index0] <- 0.2
  }
  
  return(pi_true)
  
}

#sample n times from this cdf
gensamp <- function(cdf,tgrid, n){
  unisamps <- runif(n)
  samp <- rep(0,n)
  #print(length(cdf))
  #print(length(tgrid))
  for (j in 1:n){
    index <- sum( 1*(cdf<= unisamps[j])  )
    if (index == 0){
      index = 1
    }
    samp[j] <- tgrid[ index ]
  }
  return(samp)
}

makeciplot <- function(out,sided=1,main="",baseline_factor=2.8)
{
  cs <- out$cs
  zeroindex <- which.min((out$cs - 1)^2)
  out$ests <- (out$ests-out$ests[zeroindex])
  indices <- cs >= 1
  plot(out$cs[indices], -out$ests[indices],type='l',ylim=c(0,-pnorm(1.96-sqrt(3)*1.7)+pnorm(1.96-1.7)),xlab="Sample Increase Factor",ylab="Power Gain",main=main)
  polygon(c(rev(cs[indices]),cs[indices]), -c(rev(out$cibot[indices]), out$citop[indices]),col='gray',border=NA)
  lines(out$cs[indices], -out$ests[indices],type='l')
  baseline <- pnorm(tc-baseline_factor)-pnorm(tc-sqrt(cs[indices])*baseline_factor)
  lines(out$cs[indices],baseline,col="blue")
  legend(x="topleft", legend=c('Estimate', 'All 80% cond pwr'),col=c("black", "blue"),cex=0.8,lty=1:2)
}

makeciplot_double <- function(out,out2,main="",series1="series 1",series2="series 2",baseline_factor=2.8)
{
  cs <- out$cs
  zeroindex <- which.min((out$cs - 1)^2)
  out$diffests <- (out$ests-out$ests[zeroindex])
  out2$diffests <- (out2$ests-out2$ests[zeroindex])
  indices <- cs >= 1
  plot(out$cs[indices], -out$ests[indices],type='l',ylim=c(0,0.3),xlab="Sample Increase Factor",ylab="Power Gain",main=main)
  polygon(c(rev(cs[indices]),cs[indices]), -c(rev(out$cibot[indices]), out$citop[indices]),col='gray',border=NA)
  polygon(c(rev(cs[indices]),cs[indices]), -c(rev(out2$cibot[indices]), out2$citop[indices]),col='green',border=NA)
  lines(out$cs[indices], -out$diffests[indices],type='l')
  lines(out$cs[indices], -out2$diffests[indices],type='l')
  baseline <- pnorm(tc-baseline_factor)-pnorm(tc-sqrt(cs[indices])*baseline_factor)
  lines(out$cs[indices],baseline,col="blue")

  legend(x="topleft", legend=c(series1,series2, 'All 80% cond pwr'),col=c("gray","green", "blue"),cex=0.8,lty=c(1,1,1))
}

makeciplot_double_gg <- function(out1,out2,title1, title2,main){
  zeroindex <- which.min((out1$cs - 1)^2)
  out1$diffests <- (out1$ests-out1$ests[zeroindex])
  zeroindex <- which.min((out2$cs - 1)^2)
  out2$diffests <- (out2$ests-out2$ests[zeroindex])
  fields <- data.frame(cs = c(out1$cs,out2$cs), Literature=c(rep(title1,length(out1$cs)),rep(title2,length(out2$cs))), citop = c(-out1$citop,-out2$citop), cibot = c(-out1$cibot,-out2$cibot),diffests = c(-out1$diffests,-out2$diffests) )
  ggplot(data=fields[fields$cs >=1,], aes(x=cs, y=citop,group=Literature))  +
    theme(text = element_text(size = 20)) + 
    theme(axis.text = element_text(size = 20)) + 
    geom_ribbon(aes(ymin=citop, ymax=cibot,group=Literature,fill=Literature), linetype=2, alpha=0.4) + ggtitle(main) +
    ylab("Power Gain") + xlab("Sample Size Increase Factor (c)") + 
    geom_line(aes(x=cs,y=diffests,group=Literature,color=Literature),lwd=1) + ylim(0,.3)
}

makeciplot_gg <- function(out1, main){
  zeroindex <- which.min((out1$cs - 1)^2)
  out1$diffests <- (out1$ests-out1$ests[zeroindex])
  fields <- data.frame(cs = c(out1$cs), citop = c(-out1$citop), cibot = c(-out1$cibot),diffests = c(-out1$diffests) )
  ggplot(data=fields[fields$cs >=1,], aes(x=cs, y=citop,))  +
    geom_ribbon(aes(ymin=citop, ymax=cibot), linetype=2, alpha=0.4) + ggtitle(main) +
    ylab("Power Gain") + xlab("Sample Size Increase Factor (c)") +
    geom_line(aes(x=cs,y=diffests),lwd=1)
}