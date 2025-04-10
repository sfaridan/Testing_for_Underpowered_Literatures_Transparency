
#multiple hs
rm(list = ls())
hs <- c(.1,.1,.1)
weights <- c(0.1,0.475,0.475)
weights <- weights/sum(weights)
cv <- 1.96
parm <-  2.8016177
lambda <- 1/((1-pnorm(cv-parm)+pnorm(-cv-parm))-(parm/2)*(dnorm(cv-parm) -dnorm(-cv-parm)))-1
Power <- function(ns){
  pwrs <- rep(0,length(ns))
  for (j in 1:length(ns)){
    n <- ns[j]
    pwrs[j] <- sum((1-pnorm(cv- hs*sqrt (n)) +pnorm(-cv- hs*sqrt (n))  )*weights)  
  }
  return(pwrs) 
}
U <- function(n){
  return( ( (1+lambda)*Power(n) -1 )/n )
}
maxn <- round(2*parm/min(hs))^2
ns <- 1:(maxn*10) /10
Us <- U(ns)
maximizer <- which.max(Us)
plot(Power((ns))[(maximizer-10): (maximizer+10)],U(ns[(maximizer-10): (maximizer+10)]) )
print(c(ns[maximizer] , Power(ns[maximizer]) ))


xs <- 1:600 /100
bp1 <- pnorm(cv-xs)
bp2 <- xs/2*dnorm(cv-xs)
plot(bp1,bp2)

