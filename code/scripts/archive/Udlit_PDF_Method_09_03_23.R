rm(list = ls())

library(EQL)

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
  bars2[abs(breaks)<=tc] <- bars[abs(breaks)<=tc]/theta
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


#quick unit test
theta <- 0.75
tc <- 1.96
data_pre <- rnorm(5000,mean=1,sd=2)
rs <- runif(length(data_pre))<= theta
data <- data_pre[ rs ==1 | abs(data_pre)>tc ]
limit <- 7
delta <- 0.3
breaks <- seq(tc-limit,tc+limit,delta)
histo <- histogram_fun(data,breaks)
thetahat <- histo[sum(breaks<=tc)-1]/histo[sum(breaks<=tc)]
print(thetahat)
plot(breaks[1:(length(breaks)-1)],histo,type="l")
histo_nopb <- remove_pb(breaks[2:length(breaks)],histo, thetahat,tc,delta)
plot(breaks[1:(length(breaks)-1)],histo_nopb,type="l")


#Each row of matrix is a Hermite polynomial
sigma <- 1
sigma_Y <- 2
sigma_X <- sqrt(sigma_Y^2+sigma^2)
rho <- sigma_Y^2/(sigma^2+sigma_Y^2)
J <- 12
Convmat <- matrix ( rep(0,J*length(histo)) ,nrow = J )
Deconvmat <- matrix ( rep(0,J*length(histo)) ,nrow = J )
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
  }
}
coeffs <- Convmat%*%histo_nopb
coeffs_deconv <- Deconvmat%*%histo_nopb

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
xs <- (-50:50) / 10
cva <- conv_approx_fun(xs)
plot(xs,cva,main="Convolution")
lines(xs, dnorm((xs-1)/sqrt(5))/sqrt(5),col="red")
cvd <- deconv_approx_fun(xs)
plot(xs,cvd,main="Deconvolution")
lines(xs, dnorm((xs-1)/sqrt(3))/sqrt(3),col="red")


