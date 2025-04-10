rm(list = ls())

library(EQL)



#inner product
inner_prod <- function(f1,f2,sigma_Y){
  integrand <- function(x){ return( f1(x)*f2(x)*dnorm(x/(sigma_Y))/(sigma_Y) ) }
  integral <- integrate(integrand,lower=-10, upper=10)
  return(integral$value)
}


hermite_general<- function(x,j,sigma_Y){
  ans <- 0*x
  for (l in 0:floor(j/2)){
    ans <- ans+(-1)^l*factorial(2*l)/2^l/factorial(l)*choose(j,2*l)*(x/sigma_Y)^(j-2*l)
  }
  return(ans/sqrt(factorial(j)))
}

#shifts by sqrt 2 instead of 2

#use SVD to approximate convoluting N(0,1) with N(0,1) twice (to get N(0,3))
sigma_start <- 1
fa <- function(x){return(dnorm((x-2)/sqrt(sigma_start^2))/sqrt(sigma_start^2) )}
fb <- function(x){return(dnorm((x+4)/sqrt(0.6*sigma_start^2))/sqrt(0.6*sigma_start^2) )}
f <- function(x){ return(0.5*fa(x)+0.5*fb(x))}

sigma <- 1/2
f_post <- function(x){return(dnorm((x-1)/sqrt(sigma_start^2+sigma^2))/sqrt(sigma_start^2+sigma^2))}
fa_post <-  function(x){return(dnorm((x-2)/sqrt(sigma_start^2+sigma^2))/sqrt(sigma_start^2+sigma^2) )}
fb_post <- function(x){return(dnorm((x+4)/sqrt(0.6*sigma_start^2+sigma^2))/sqrt(0.6*sigma_start^2+sigma^2) )}
f_post <- function(x){ return(0.5*fa_post(x)+0.5*fb_post(x))}

sigma_Y <- 1
sigma_X <- sqrt(sigma_Y^2+sigma^2)
rho <- sigma^2/(sigma^2+sigma_Y^2)
rho <- sigma_Y^2/(sigma^2+sigma_Y^2)
numsvs <- 16
coeffs <- rep(NA,numsvs)
for (j in 1:numsvs){
  lambda <- rho^((j-1)/2) #*(1/4)^((j-1)/2)
  hj <- function(x){return(hermite_general(x,j-1,sigma_X))}
  print(inner_prod(hj,hj,1))
  coeffs[j] <- lambda * inner_prod(f,hj,sigma_X) #/inner_prod(singular_function,singular_function,sigma_Y)
}

coeffs_deconv <- rep(NA,numsvs)
for (j in 1:numsvs){
  lambda <- rho^((j-1)/2) #*(1/4)^((j-1)/2)
  hj <- function(x){return(hermite_general(x,j-1,sigma_Y))}
  print(inner_prod(hj,hj,1))
  coeffs_deconv[j] <- (1/lambda) * inner_prod(f_post,hj,sigma_Y) #/inner_prod(singular_function,singular_function,sigma_Y)
}

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
decva <- deconv_approx_fun(xs)

plot(xs,decva)
lines(xs,f(xs),type="l",col="red")

plot(xs,cva)
lines(xs,f_post(xs),type="l",col="red")

