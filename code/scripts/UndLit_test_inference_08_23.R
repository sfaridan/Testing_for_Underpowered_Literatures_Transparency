
rm(list = ls())
library(Rsolnp)
library(matlib)

empcdf_bygrid<- function(Tscore_data,tgrid){
  
  ecdf <- 0*tgrid
  for (t in 1:length(ecdf)){
    ecdf[t] = mean(Tscore_data<= tgrid[t])
  }
  return(ecdf)
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

#returns gradient of D(Bnpi,Theta)
Gradient <- function(Bn,pis,Theta,tgrid,tc){
  G <- matrix(rep(0,dim(Bn)[1]*(dim(Bn)[2]+1)),nrow = dim(Bn)[1])
  index <- sum(tgrid<=tc)
  Bnpi <- Bn%*%pis
  td <- tgrid<=tc
  tu <- tgrid>tc
  denom <- 1-(1-Theta)*(Bn[index,]%*%pis)
  numer <- Theta*Bnpi*td + Theta*Bnpi[index]*tu + (Bnpi-Bnpi[index])*tu
  
  for (i in 1:(dim(Bn)[1]) ){
    for (j in 1:(dim(Bn)[2]) ){
      G[i,j] <-  (Theta*Bn[i,j]*td[i] +Theta*Bn[index,j]*tu[i]+(Bn[i,j]-Bn[index,j])*tu[i] )/denom + (1-Theta)*Bn[index,j]*(numer[i])/denom^2
      }
    G[i,dim(G)[2]] <- (Bnpi[i]*td[i]+Bnpi[index]*tu[i])/denom-Bnpi[index]*numer[i]/denom^2
  }
  return(G)
}

getOmega <- function(n,Fhat){
  Omega <- matrix( rep(0,length(Fhat)^2),nrow=length(Fhat))
  for (i in 1:length(Fhat)){
    for(j in i:length(Fhat)){
      Omega[i,j] <- Fhat[i]*abs(1-Fhat[j])/n
      Omega[j,i] <- Omega[i,j]
    }
  }
  return(Omega)
}

pub_bias_matrix <- function(tgrid, mat,tc,Theta){
  
  mat_withpb <- 0*mat
  for (i in 1:(dim(mat)[2])){
    mat_withpb[,i] <-  pub_bias(tgrid, mat[,i],tc,Theta)
  }
  
  return(mat_withpb)
}

pub_bias <- function(tgrid, cdf,tc,Theta){
  ts <- length(tgrid)
  index <- sum(tgrid<=tc)
  denom <- 1-(1-Theta)*cdf[index]
  cdf_withpb <-  (  Theta*cdf*(tgrid<= tc)+Theta*cdf[index]*(tgrid> tc)+(cdf-cdf[index])*(tgrid> tc) )/ denom
  return(cdf_withpb)
}

estimator <- function(ecdf,Bn,tgrid,tc){
  obj<- function(X){
    pihat <- X[1:(length(X)-1)]
    thetahat <- X[length(X)]
    resids <- pub_bias(tgrid, Bn%*%pihat, tc, thetahat) - ecdf
    #resids <- Bn%*%pihat - ecdf
    return(100000*t(resids)%*%resids )
  }

  hl <- dim(Bn)[2]
  X0 <- c(rep(1/hl,hl),0.9)
  sols <- solnp(X0, obj,LB =c(rep(-100,hl),0.01),control = list(trace = 0))
  ests <- list()
  ests$pihat <- sols$pars[1:(dim(Bn)[2])]
  ests$thetahat <- sols$pars[length(sols$pars)]
  ests$fit <-obj(sols$pars)
  ests$hessian <- sols$hessian
  return(ests)
}


estimator_grad <- function(ecdf,Bn,tgrid,tc,G){
  obj<- function(X){
    pihat <- X[1:(length(X)-1)]
    thetahat <- X[length(X)]
    resids <- t(G)%*%(pub_bias(tgrid, Bn%*%pihat, tc, thetahat) - ecdf)
    #resids <- Bn%*%pihat - ecdf
    return(1000000*t(resids)%*%resids )
  }
  
  hl <- dim(Bn)[2]
  X0 <- c(rep(1/hl,hl),0.9)
  sols <- solnp(X0, obj,LB =c(rep(-100,hl),0.01),control = list(trace = 0))
  ests <- list()
  ests$pihat <- sols$pars[1:(dim(Bn)[2])]
  ests$thetahat <- sols$pars[length(sols$pars)]
  ests$fit <-obj(sols$pars)
  ests$hessian <- sols$hessian
  print(sols$convergence)
  return(ests)
}

estimator_long <- function(ecdf,Bn,tgrid,tc){
  obj<- function(thetahat){

    Bn_withpb <- pub_bias_matrix(tgrid, Bn, tc, thetahat)
    #coeffs <- inv(t(Bn_withpb)%*%Bn_withpb)%*%t(Bn_withpb)%*%ecdf
    coeffs <- inv(t(Bn_withpb)%*%Bn_withpb)%*%t(Bn_withpb)%*%ecdf
    resids <- Bn_withpb%*%coeffs-ecdf
    
    return(100000*t(resids)%*%resids )
  }
  
  X0 <- c(0.9)
  sols <- solnp(X0, obj,LB =c(0.01),control = list(trace = 0))
  ests <- list()
  ests$thetahat <- sols$pars[length(sols$pars)]
  
  Bn_withpb <- pub_bias_matrix(tgrid, Bn, tc,  ests$thetahat )
  coeffs <- inv(t(Bn_withpb)%*%Bn_withpb)%*%t(Bn_withpb)%*%ecdf
  ests$pihat <- inv(t(Bn)%*%Bn)%*%t(Bn)%*%pub_bias(tgrid,Bn_withpb%*%coeffs,tc,1/ests$thetahat)
  
  ests$fit <-obj(sols$pars)
  ests$hessian <- sols$hessian
  return(ests)
}



var_estimator <- function(Bn,pi_hat,Theta_hat,tgrid,tc,n,F_hat,vvc){
  tcindex <- sum(tgrid<=tc)
  G <- Gradient(Bn,pi_hat,Theta_hat,tgrid,tc)
  GpG <- t(G)%*%G
  regmat <- inv(GpG)%*%t(G)
  varvec <- vvc%*%regmat
  Omega <- getOmega(n,F_hat)
  varhat <- varvec%*%Omega%*%t(varvec)
  return(varhat)
}

tc <- 1.65
ngrid <- 15
tol <- 0.005
tgrid <- -100:100/20
tgrid <- sort(c(tgrid, tc+0.001, tc-0.001))
tcindex <- sum(tgrid<=tc)
hgrid <- -20:20/10
Bn_pre    <- cdf_mat(tgrid,hgrid,c=1,sided=1)
Bn_svd <- svd(Bn_pre)
norm(Bn_svd$u%*%diag(Bn_svd$d)%*%t(Bn_svd$v)-Bn_pre)
keeps <- abs(Bn_svd$d)>tol
Lambda <- diag(Bn_svd$d[keeps])
Bn_V <- Bn_svd$v[,keeps]
Bn <- Bn_svd$u[,keeps]%*%Lambda

hmass <- 1
pi_orig <- rep(0,length(hgrid))
pi_orig[sum(hgrid<=hmass)] <- 1
pi_reg <- inv(t(Bn_svd$v)%*%(Bn_svd$v))%*%t(Bn_svd$v)%*%pi_orig
pi_true <- pi_reg[keeps]
norm((Bn_svd$v)%*%pi_reg-pi_orig)
norm(Bn_pre%*%pi_orig - Bn_pre%*%Bn_svd$v%*%pi_reg)
norm(Bn_pre%*%pi_orig - Bn_svd$u%*%diag(Bn_svd$d)%*%pi_reg)
norm(Bn_pre%*%pi_orig - Bn_svd$u%*%diag(Bn_svd$d*(abs(Bn_svd$d)>tol))%*%pi_reg)
norm(Bn_pre%*%pi_orig - Bn_svd$u[,keeps]%*%diag(Bn_svd$d[keeps])%*%(pi_true))
norm(Bn_pre%*%pi_orig - Bn%*%(pi_true))



Theta_true <- 0.5
F_true <-  pub_bias(tgrid, Bn%*%pi_true, tc, Theta_true)
Vc     <- diag(length(tgrid)+1)
plot(tgrid,Bn%*%pi_true,type="l")
lines(tgrid,F_true,col="blue")
norm(F_true - pub_bias(tgrid, Bn_pre%*%pi_orig, tc, Theta_true))

Tscore_data <- c(rnorm(n,hmass,1)) 
rand <- runif(length(Tscore_data))
Tscore_data <- Tscore_data[(Tscore_data) >= rep(tc,length(Tscore_data)) | rand <= Theta_true ]
#Tscore_data <- (Tscore_data)
ecdf <- empcdf_bygrid(Tscore_data,tgrid)
plot(tgrid,F_true)
lines(tgrid, ecdf)


#use truth as estimate just to see how reasonable asymptotic variance is
pi_hat <- pi_true
F_hat   <- F_true
Theta_hat <- Theta_true
n      <- 500000
vvc <- c(Bn[tcindex,],0)

#simulations
nsims <- 1000
ests <- rep(0,nsims)
thetahats <- rep(0,nsims)
gradients <- rep(0,nsims)
varests <- rep(0,nsims)
gradients_truegrad <- rep(0,nsims)
gradients_Ftrue <- rep(0,nsims)
diffs <- rep(0,nsims)
covers <- rep(0,nsims)
Gtrue <- Gradient(Bn,pi_true,Theta_true,tgrid,tc)
Omega_true <- Omega(n,F_true)
regmat_true <- vvc%*%inv(t(Gtrue)%*%Gtrue)%*%t(Gtrue)
varhat_trueF <- var_estimator(Bn,pi_true,Theta_true,tgrid,tc,n,F_true,c(vvc))
print(c(sqrt(varhat_trueF) ))
for (sim in 1:nsims){
  Tscore_data <- c(rnorm(n,hmass,1)) 
  rand <- runif(length(Tscore_data))
  Tscore_data <- Tscore_data[(Tscore_data) >= rep(tc,length(Tscore_data)) | rand <= Theta_true ]
  #Tscore_data <- (Tscore_data)
  ecdf <- empcdf_bygrid(Tscore_data,tgrid)
  
  #estimate <- estimator(ecdf,Bn,tgrid,tc)
  estimate <- estimator_grad(ecdf,Bn,tgrid,tc, Gtrue)
  G <- Gradient(Bn,estimate$pihat,estimate$thetahat,tgrid,tc)
  
  
  gradients[sim]<- sqrt(n)*max(abs(t(G)%*%(pub_bias(tgrid,Bn%*%estimate$pihat,tc,estimate$thetahat)-ecdf) )) #/max(abs((pub_bias(tgrid,Bn%*%estimate$pihat,tc,estimate$thetahat)-ecdf))) 
  gradients_Ftrue[sim]<- sqrt(n)*max(abs(t(G)%*%(pub_bias(tgrid,Bn%*%pi_true,tc,Theta_true)-ecdf) )) #/max(abs((pub_bias(tgrid,Bn%*%estimate$pihat,tc,estimate$thetahat)-ecdf)))
  gradients_truegrad[sim]<- sqrt(n)*max(abs(t(Gtrue)%*%(pub_bias(tgrid,Bn%*%estimate$pihat,tc,estimate$thetahat)-ecdf) )) #/max(abs((pub_bias(tgrid,Bn%*%estimate$pihat,tc,estimate$thetahat)-ecdf)))
  
  
  
  ests[sim] <- vvc%*%c(estimate$pihat,estimate$thetahat)
  thetahats[sim] <- estimate$thetahat
  
  diffest <- regmat_true%*%(F_true-ecdf)
  diff <- vvc%*%c(pi_true,Theta_true)-ests[sim]
  diffs[sim] <- abs((diffest-diff)/diffest)
  #ih <- vvc%*%inv(estimate$hessian)
  varests[sim] <- var_estimator(Bn,estimate$pihat,estimate$thetahat,tgrid,tc,n,ecdf,c(vvc))
  
  covers[sim] <-  ( vvc%*%c(pi_true,Theta_true) <= ests[sim] + 1.96*sqrt( varests[sim]) & vvc%*%c(pi_true,Theta_true) >= ests[sim] - 1.96*sqrt( varests[sim])   )
  
  print("")
  print(c(sim, nsims))
  print("gradients: ")
  print(c(mean(gradients[1:sim]),mean(gradients_truegrad[1:sim]), sum(((pub_bias(tgrid,Bn%*%estimate$pihat,tc,estimate$thetahat)-ecdf))) ))
  print("differences")
  print(c(mean(diffs[1:sim])))
  print(c(mean(thetahats[1:sim]), Theta_true, sd(thetahats[1:sim])))
  print(c(mean(ests[1:sim]), vvc%*%c(pi_true,Theta_true)))
  print(sqrt(n)*c(sd(ests[1:sim]), sqrt(varhat_trueF),mean(sqrt(varests[1:sim]))))
  print(mean(covers[1:sim]))
}

#check if we have the actual gradient
th <- .25
Gtrue <- Gradient(Bn,pi_true,th,tgrid,tc)
j <- 1
eps <- abs(pi_true[1]-estimate$pihat[1])
delta <- 0*pi_true
delta[j] <- eps
F1 <-pub_bias(tgrid,Bn%*%pi_true,tc,th)
F2 <-pub_bias(tgrid,Bn%*%(pi_true+delta),tc,th)
cbind((F2-F1)/eps,Gtrue[,j])
print(sum(abs(  (F2-F1)/eps-Gtrue[,j] )))
print(sum(abs(  (F2[140]-F1[140])/eps-Gtrue[140,j] )))
t(  (F2-F1)/eps-Gtrue[,j] )%*%(F_true-ecdf)

F1th <-pub_bias(tgrid,Bn%*%pi_true,tc,th)
F2th <-pub_bias(tgrid,Bn%*%(pi_true),tc,th+eps)
cbind((F2th-F1th)/eps,Gtrue[,dim(G)[2]])
print(sum(abs( (F2th-F1th)/eps-Gtrue[,dim(G)[2]] )))

