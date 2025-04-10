rm(list = ls())

library(EQL)

pdf_mat <-  function(tgrid,hgrid){
  
  argst = matrix(rep(tgrid,length(hgrid)),ncol=length(hgrid))
  argsh = t(matrix(rep(hgrid,length(tgrid)),ncol=length(tgrid)))
  
  V = dnorm(argst - argsh)
  
  return(V)
}

#returns a function with a histogram made from the data 
# bar breaks are delta apart, center on tc and go from tc-limit to tc+limit
histogram_fun <- function(data,breaks){
  bars <- 0*rep(0,length(breaks)-1)
  for (j in 1:(length(breaks)-1)){
    bars[j] <- mean(data < breaks[j+1] & data >= breaks[j])/delta
  }
  return(bars)
}

remove_pb <- function(breaks,bars,theta,tc,delta){
  bars2 <- bars
  includes <- abs(breaks)<=tc 
  bars2[includes] <- bars[includes]/theta
  bars2 <- bars2 / (sum(bars2))/delta
  return(bars2)
}

#inner product
inner_prod <- function(f1,f2,sigma_Y){
  integrand <- function(x){ return( f1(x)*f2(x)*dnorm(x/(sigma_Y))/(sigma_Y) ) }
  integral <- integrate(integrand,lower=-5*sigma_Y, upper=5*sigma_Y,subdivisions=2000)
  return(integral$value)
}


hermite_general<- function(x,j,sigma_Y){
  ans <- 0*x
  for (l in 0:floor(j/2)){
    ans <- ans+(-1)^l*factorial(2*l)/2^l/factorial(l)*choose(j,2*l)*(x/sigma_Y)^(j-2*l)
  }
  return(ans/sqrt(factorial(j)))
}

#preliminaries
theta <- 1
tc <- 1.96
n <- 78000
limit <- 4
delta <- (2*(tc+limit))/60
breaks <- seq(-tc-limit,tc+limit,delta)
tops <- breaks[1:(length(breaks)-1)]
hgrid <- seq(-tc-2,tc+2,0.5)
tgrid_l <- tops[abs(tops)<tc-2*delta ]
tgrid_u <- tops[abs(tops)>=tc+2*delta]
Bl <- pdf_mat(tgrid_l, hgrid)
Bu <- pdf_mat(tgrid_u, hgrid)
#Reg_L <- inv(t(Bl)%*%Bl +.01*diag(dim(Bl)[2]))%*%t(Bl)
#Reg_U <- inv(t(Bu)%*%Bu +.01*diag(dim(Bu)[2]))%*%t(Bu)


nsims <- 200
hats <- rep(0,nsims)
changehat <- rep(0,nsims)
thetahats <- rep(0,nsims)
for (sim in 1:nsims){

data_pre <- rnorm(n,mean=0,sd=2)
rs <- runif(length(data_pre))<= theta
data <- data_pre[ rs ==1 | abs(data_pre)>tc ]
histo <- histogram_fun(data,breaks)
histo <- histo /sum(histo) * delta

thetahat_old <- mean(  abs(data)< 1.96 & abs(data)>(1.96-.2/1.77^3)  )/ mean(  abs(data)>= 1.96 & abs(data)< (1.96+.2/1.77^3)   )  #(histo[sum(breaks<=tc)]+histo[sum(breaks<=-tc)])/(histo[sum(breaks<=tc)+1]+histo[sum(breaks<=-tc)+1])

#new thetahat
#pihat_l <- Reg_L%*%histo[abs(tops)<tc-2*delta ]
#pihat_u <- Reg_U%*%histo[abs(tops)>=tc+2*delta ]
#v <- pdf_mat(tc,hgrid)+pdf_mat(-tc,hgrid)
#pl <- sum(v*t(pihat_l))
#pu <- sum(v*t(pihat_u))
#thetahat <- pl/pu

#plot(breaks[abs(tops)<tc-2*delta], histo[abs(tops)<tc-2*delta] ,type="l")
#points(breaks[abs(tops)<tc-2*delta],Bl%*%pihat_l )
#plot(tops[abs(tops)>=tc+2*delta], histo[abs(tops)>=tc+2*delta] ,type="l")
#points(tops[abs(tops)>=tc+2*delta],Bu%*%pihat_u )
#0.001979034
#[1] "sd est:0.0116341670725002 sd theta: 0.0448277060035314"

thetahats[sim] <- theta #thetahat_old #theta +rnorm(1)*0.28+.1
#thetahat <- min(c(thetahat,1.01))
plot(breaks[1:(length(breaks)-1)],histo,type="l")
histo_nopb <- remove_pb(breaks[2:length(breaks)],histo, thetahats[sim],tc,delta)
plot(breaks[1:(length(breaks)-1)],histo_nopb,type="l")


#Each row of matrix is a Hermite polynomial
sigma <- 1/2
sigma_Y <- 3
sigma_X <- sqrt(sigma_Y^2+sigma^2)
rho <- sigma_Y^2/(sigma^2+sigma_Y^2)
J <- 20
J <- ceiling(2*log(1/n^(2/3))/log(sigma_Y^2/(1+sigma_Y^2)))
Convmat <- matrix ( rep(0,J*length(histo)) ,nrow = J )
Deconvmat <- matrix ( rep(0,J*length(histo)) ,nrow = J )
Origmat <- matrix ( rep(0,J*length(histo)) ,nrow = J )
for (j in 1:J){
  lambda <- rho^((j-1)/2) 
  hj <- function(x){return(hermite_general(x,j-1,sigma_X))}
  hde <- function(x){return(hermite_general(x,j-1,sigma_Y))}
  integrand <- function(x){ return( hj(x)*dnorm(x/(sigma_X))/(sigma_X) ) }
  integrand_de <- function(x){ return( hde(x)*dnorm(x/(sigma_Y))/(sigma_Y) ) }
  for (br in 1:length(histo)){
    lb <- breaks[br]
    ub <- breaks[br+1]
    Convmat[j,br] <- lambda*integrate(integrand,lower=lb, upper=ub,subdivisions=2000)$value
    Deconvmat[j,br] <- (1/lambda)*integrate(integrand_de,lower=lb, upper=ub,subdivisions=2000)$value
    Origmat[j,br] <- integrate(integrand_de,lower=lb, upper=ub,subdivisions=2000)$value
  }
}
coeffs <- Convmat%*%histo_nopb
coeffs_deconv <- Deconvmat%*%histo_nopb
coeffs_orig <- Origmat%*%histo_nopb

conv_approx_fun <- function(xs){
  ans <- 0*xs
  for (j in 1:length(coeffs)){
    hj <- function(x){return(hermite_general(x,j-1,sigma_Y))}
    for (k in 1:length(xs)){
      ans[k] <- ans[k] + coeffs[j]*hj(xs[k])
    }
  }
  return(ans)
}
deconv_approx_fun <- function(xs){
  ans <- 0*xs
  for (j in 1:length(coeffs)){
    hj <- function(x){return(hermite_general(x,j-1,sigma_X))}
    for (k in 1:length(xs)){
      ans[k] <- ans[k] + coeffs_deconv[j]*hj(xs[k])
    }
  }
  return(ans)
}
orig_approx_fun <- function(xs){
  ans <- 0*xs
  for (j in 1:length(coeffs)){
    hj <- function(x){return(hermite_general(x,j-1,sigma_X))}
    for (k in 1:length(xs)){
      ans[k] <- ans[k] + coeffs_orig[j]*hj(xs[k])
    }
  }
  return(ans)
}
xs <- (-50:50) / 10
cva <- conv_approx_fun(xs)
#plot(xs,cva,main="Convolution with pb removed")
#lines(xs, dnorm((xs-0)/sqrt(2^2+1/2^2))/sqrt(2^2+1/2^2),col="red")
cvd <- deconv_approx_fun(xs)
#plot(xs, dnorm((xs-0)/sqrt(2^2-1/2^2))/sqrt(2^2-1/2^2),col="red")
#lines(xs,cvd,main="Deconvolution with pb removed")

integrand_orig <- function(x){dnorm((x)/sqrt(2^2))/sqrt(2^2)}
integrand_true <- function(x){dnorm((x)/sqrt(2^2-1/2^2))/sqrt(2^2-1/2^2)}
integrand_hat <- function(x){deconv_approx_fun(x)}
integrand_hat_orig <- function(x){orig_approx_fun(x)}
Forig <- integrate(integrand_orig, -tc,tc)$value
Ftrue <- integrate(integrand_true, -tc/sqrt(2),tc/sqrt(2))$value 
Fhat <- integrate(integrand_hat, -tc/sqrt(2),tc/sqrt(2))$value 
Fhatorig <- integrate(integrand_hat_orig, -tc,tc)$value
changehat[sim] <- Fhat-Fhatorig
hats[sim] <- Fhat 
print("hats::")
print(c(Ftrue,Fhat))
print(c(Forig,Fhatorig))
print(c(Ftrue-Forig, mean(changehat[1:sim])))
print(c(sd(hats[1:sim]),sd(changehat[1:sim])))
print(paste0("sd est:" ,sd(hats[1:sim]), " sd theta: ", sd(thetahats[1:sim])))
print(paste0("true: " ,Ftrue, " mean: ", mean(hats[1:sim])))
}
hist(hats[1:sim-1])
print(sd(hats[1:sim-1]))
print(c(Ftrue,mean(hats[1:sim-1])))

