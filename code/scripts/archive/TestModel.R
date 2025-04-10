
#multiple hs
rm(list = ls())
hs <- c(.11,.1,.1)
weights <- c(0.1,0.475,0.475)
weights <- weights/sum(weights)
cv <- 1.96
beta <- 0.2
parm <- cv-qnorm(beta)
lambda <- 1/(1-beta-(parm/2)*dnorm(cv-parm))-1
Power <- function(ns){
  pwrs <- rep(0,length(ns))
  for (j in 1:length(ns)){
    n <- ns[j]
    pwrs[j] <- sum((1-pnorm(cv- hs*sqrt (n))  )*weights)  
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


xs <- -500:600 /100
bp1 <- pnorm(cv-xs)
bp2 <- xs/2*dnorm(cv-xs)
plot(bp1,bp2)

plot(bp1[1:200],bp2[1:200],type="l")

cv <- 1.65
pnorm(0.5*(cv-sqrt(cv^2+4)))
bs <- 1:999 / 1000
vals<-(cv-qnorm(bs))/2*dnorm(qnorm(bs))-cv/4*dnorm(cv/2)
print(cbind(bs,sign(vals))  )

cvs <- 0:100/10
vals<-(cvs-qnorm(0.3))/2*dnorm(qnorm(0.3))-cvs/4*dnorm(cvs/2)
print(cbind(cvs,sign(vals))  )