
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
  reg <- inv(t(Bn)%*%Bn)%*%t(Bn)
  proj <- Bn%*%reg
  
  obj<- function(thetahat){
    ecdf_de_pb <- pub_bias(tgrid, ecdf, tc, 1/thetahat)
    resids <- lm.fit(Bn,ecdf_de_pb)$residuals
    #resids <- proj%*%ecdf_de_pb-ecdf_de_pb
    #resids <- Bn%*%coeffs-ecdf_de_pb
    return(t(resids)%*%resids )
    #G <- Gmat(Bn,thetahat, ecdf, tgrid, tc)
    #resids <- t(G)%*%(proj%*%ecdf_de_pb-ecdf_de_pb)
    #max(abs(resids))
  }
  
  hl <- dim(Bn)[2]
  X0 <- c(0.9)
  sols <- solnp(X0, obj,LB =c(0.01),UB=2,control = list(trace = 0))
  ests <- list()
  ests$thetahat <- sols$pars[length(sols$pars)]
  ecdf_de_pb <- pub_bias(tgrid, ecdf, tc, 1/ests$thetahat)
  ests$pihat <- lm.fit(Bn,ecdf_de_pb)$coefficients
  ests$fit <-obj(sols$pars)
  ests$hessian <- sols$hessian
  return(ests)
}





Gmat <- function(Bn,theta, Fhat, tgrid, tc){
  theta <- 1/theta #ALERT
  
  theta_column <-rep(0,dim(Bn)[1])
  td <- tgrid <= tc
  tu <- tgrid > tc
  Ftc <- Fhat[sum(tgrid<=tc)]
  denom <- 1-(1-theta)*Ftc
  
  for (i in 1:(dim(Bn)[1])){
    term1 <- (Fhat[i]*td[i]+Ftc*tu[i])/denom
    term2 <- -Ftc*(theta*Fhat[i]*td[i] + theta*Ftc*tu[i]+ tu[i]*(Fhat[i]-Ftc))/(denom^2)
    theta_column[i]  <- term1+term2
  }
  G <- cbind(Bn, theta_column)
  return(G)
}


tc <- 1.65
ngrid <- 15
tol <- 0.0001
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

hmass <- 0.5
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



Theta_true <- 0.75
F_true <-  pub_bias(tgrid, Bn%*%pi_true, tc, Theta_true)
Vc     <- diag(length(tgrid)+1)
plot(tgrid,Bn%*%pi_true,type="l")
lines(tgrid,F_true,col="blue")
#check that dim reduction hasn't changed anything
norm(F_true - pub_bias(tgrid, Bn_pre%*%pi_orig, tc, Theta_true))


n      <- 78
vvc <- c(Bn[tcindex,],0)

#population
Tscore_pop<- c(rnorm(1000000,hmass,1)) 
rand <- runif(length(Tscore_pop))
Tscore_pop <- Tscore_pop[(Tscore_pop) >= rep(tc,length(Tscore_pop)) | rand <= Theta_true ]


#check that population gradient is zero
Gpop <- Gmat(Bn,Theta_true , F_true, tgrid, tc)
gradpop <- t(Gpop)%*%(  Bn%*%pi_true-pub_bias(tgrid,F_true,tc,1/Theta_true) )
print(t(gradpop))

#compute asymptotic variance
Omega_pop <- getOmega(n,F_true)
varestmat_pop <- inv(t(Gpop)%*%Gpop)%*%t(Gpop)
Apop <- 0*Omega_pop
inth <- 1/Theta_true
denom <- (1-(1-inth)*F_true[sum(tgrid<=tc)])
for (ai in 1:dim(Omega_pop)[1]){
  Apop[ai,ai] <- (inth*(tgrid[ai]<=tc) +(tgrid[ai]>tc))/denom
  Apop[ai, sum(tgrid<=tc)] <- Apop[ai, sum(tgrid<=tc)] + (inth-1)*(tgrid[ai]>tc)/denom + (1-inth)*(inth*F_true[ai]*(tgrid[ai]<=tc) +inth*F_true[sum(tgrid<=tc)]*(tgrid[ai]>tc)+(tgrid[ai]>tc)*(F_true[ai]-F_true[sum(tgrid<=tc)]) )/denom^2    
}
avar <- varestmat_pop[1,]%*%Apop%*%Omega_pop%*%t(Apop)%*%varestmat_pop[1,]
print(avar)

#simulations
nsims <- 2000
ests <- rep(0,nsims)
thetahats <- rep(0,nsims)
gradients <- rep(0,nsims)
varests <- rep(0,nsims)
gradients_truegrad <- rep(0,nsims)
gradients_Gpop <- rep(0,nsims)
diffs <- rep(0,nsims)
derrs <- rep(0,nsims)
derrpop <- rep(0,nsims)
cover <- rep(0,nsims)
tester <- rep(0,nsims)
tester2 <- rep(0,nsims)
Omega_true <- Omega(n,F_true)
for (sim in 1:nsims){
  
  Tscore_data <- sample(Tscore_pop, n,replace=TRUE )
  ecdf <- empcdf_bygrid(Tscore_data,tgrid)
  
  #estimate <- estimator(ecdf,Bn,tgrid,tc)
  estimate <- estimator(ecdf,Bn,tgrid,tc)
  
  #ests[sim] <- vvc%*%c(estimate$pihat,estimate$thetahat)
  ests[sim] <- estimate$pihat[1]
  thetahats[sim] <- estimate$thetahat
 
  
  G <- Gmat(Bn,thetahats[sim] , ecdf, tgrid, tc)
  gradients[sim] <- max(abs( t(G)%*%(Bn%*%estimate$pihat - pub_bias(tgrid,ecdf,tc,1/estimate$thetahat) )))
  gradients_Gpop[sim] <- max(abs( t(Gpop)%*%(Bn%*%estimate$pihat - pub_bias(tgrid,ecdf,tc,1/estimate$thetahat) )))
 
  
  #estimate the variance
  Omega <- getOmega(n,ecdf)
  varestmat <- inv(t(G)%*%G)%*%t(G)
  Amat <- 0*Omega
  inth <- 1/thetahats[sim]
  denom <- (1-(1-inth)*ecdf[sum(tgrid<=tc)])
  for (ai in 1:dim(Omega)[1]){
    Amat[ai,ai] <- (inth*(tgrid[ai]<=tc) +(tgrid[ai]>tc))/denom
    Amat[ai, sum(tgrid<=tc)] <- Amat[ai, sum(tgrid<=tc)] + (inth-1)*(tgrid[ai]>tc)/denom + (1-inth)*(inth*ecdf[ai]*(tgrid[ai]<=tc) +inth*ecdf[sum(tgrid<=tc)]*(tgrid[ai]>tc)+(tgrid[ai]>tc)*(ecdf[ai]-ecdf[sum(tgrid<=tc)]) )/denom^2    
  }
  
  #This exercise proves that the problem is with Omega
  varests[sim] <- (varestmat%*%Amat%*%Omega%*%t(Amat)%*%t(varestmat))[1,1]
  #derrs[sim] <- (inv(t(G)%*%G)%*%t(G)%*%Amat%*%(F_true-ecdf))[1]
  #derrpop[sim] <- (inv(t(Gpop)%*%Gpop)%*%t(Gpop)%*%Apop%*%(F_true-ecdf))[1]
  
  #tester[sim] <- sum(Tscore_data <= tgrid[120])/length(Tscore_data)
  #tester2[sim] <- sum(Tscore_data <= tgrid[120])/elength
  
  cover[sim] <- (ests[sim] <=pi_true[1] + 1.96*sqrt(varests[sim])) & (ests[sim] >=pi_true[1] - 1.96*sqrt(varests[sim]))
  
  print("")
  print(c(sim, nsims))
  print(c(mean(thetahats[1:sim]), Theta_true, sd(thetahats[1:sim])))
  print(c(mean(ests[1:sim]), pi_true[1],sd(ests[1:sim])))
  print((n)*c((mean(varests[1:sim])),avar,sd(ests[1:sim])^2))
  print(mean(cover[1:sim]))
#  print(c(mean(gradients[1:sim]),mean(gradients_Gpop[1:sim]) ))
  #print(c( length(Tscore_data)/elength))
  #print(c(sd(tester[1:(sim-1)])^2,(F_true[120]*(1-F_true[120]))/elength,(ecdf[120]*(1-ecdf[120]))/length(Tscore_data)))
  #print(c(sd(tester[1:sim])^2,sd(tester2[1:sim])^2,(Bn%*%pi_true)[120]*Theta_true*(1-(Bn%*%pi_true)[120]*Theta_true)/( length(Tscore_data))))
  
  #print(sqrt(n)*c(sd(ests[1:sim]), sqrt(varhat_trueF),mean(sqrt(varests[1:sim]))))
  #print(mean(covers[1:sim]))
}
  
#Check taylor approximation of the de-pb functional
sqrt(n)*cbind(Amat%*%(ecdf-F_true),  pub_bias(tgrid, ecdf, tc, 1/Theta_true) -pub_bias(tgrid, F_true, tc, 1/Theta_true)  )

#Check taylor approximation of the distance function
diffgrad <-((Bn%*%pi_true - pub_bias(tgrid,ecdf,tc,1/Theta_true) ))-((Bn%*%estimate$pihat - pub_bias(tgrid,ecdf,tc,1/estimate$thetahat) ))
diff2 <- -G%*%((c(estimate$pihat,estimate$thetahat)- c(pi_true,Theta_true)))
diff3 <- -Gpop%*%(c(estimate$pihat,estimate$thetahat)- c(pi_true,Theta_true))
cbind(diffgrad,diff2,diff3)

#establish that estimation error can be written as regression on g
err <- c(estimate$pihat,estimate$thetahat)- c(pi_true,Theta_true)
derr <- -inv(t(G)%*%G)%*%t(G)%*%(Bn%*%pi_true - pub_bias(tgrid,ecdf,tc,1/Theta_true) )
derr2 <- -inv(t(G)%*%G)%*%t(G)%*%(pub_bias(tgrid,F_true,tc,1/Theta_true) - pub_bias(tgrid,ecdf,tc,1/Theta_true) )
derr3 <- -inv(t(Gpop)%*%Gpop)%*%t(Gpop)%*%Apop%*%(F_true-ecdf)
derr4 <- -inv(t(G)%*%G)%*%t(G)%*%Amat%*%(F_true-ecdf)
sqrt(n)*cbind(err, derr,derr2,derr3,derr4) #wait this is working??

#Experiment: finding theta by hand
n <- 1000
Theta_true <- 0.5
Tscore_data <- c(rnorm(n,hmass,1)) 
rand <- runif(length(Tscore_data))
Tscore_data <- Tscore_data[(Tscore_data) >= rep(tc,length(Tscore_data)) | rand <= Theta_true ]
ecdf <- empcdf_bygrid(Tscore_data,tgrid)
thetagrid <- 100:990 / 1000
dists <- rep(0,length(thetagrid))
grads <- rep(0,length(thetagrid))
for (i in 1:length(thetagrid)){
  ecdf_de_pb <- pub_bias(tgrid, ecdf, tc, 1/thetagrid[i])
  resids <- lm.fit(Bn,ecdf_de_pb)$residuals
  dists[i] <- t(resids)%*%resids
  G <- Gmat(Bn,thetagrid[i] , ecdf, tgrid, tc)
  Gvec <- t(G)%*%resids
  grads[i] <- sum(abs(Gvec))
}
plot(thetagrid,(dists-min(dists))/mean(dists),type="l")
lines(thetagrid, grads/mean(grads),col="blue")
abline(v=Theta_true,col="red")
