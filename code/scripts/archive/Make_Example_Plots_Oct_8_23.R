ns <- (10:30) *10
sims <- 10000
pwrs1 <- 0*ns
pwrs2 <- 0*ns
library(latex2exp )


hs1 <- rnorm(sims,1.7/10,0.1 )
hs2 <- c(rnorm(sims/2,.2/10,0.1 ), rnorm(sims/2,20/10,0.1))
ts  <- rnorm(sims)

for (nn in 1:length(ns)){
  n <- ns[nn]
  pwrs1[nn] <- mean((sqrt(n)*(hs1)+ts) >1.96)
  pwrs2[nn] <- mean((sqrt(n)*(hs2)+ts) >1.96)
}


setwd("C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/output/figures" )

pdf(file = "example_two_literatures.pdf",   # The directory you want to save the file in
    width = 6*2.2 , # The width of the plot in inches *#4*1.6
    height = 6)
plot(ns/100, pwrs1,type="line",cex.lab=1.3,lwd=3, ylab=TeX(r'(Unconditional Power)'), xlab=TeX(r'(Sample Size Increase: $c^2$)'))
lines(ns/100,pwrs2,lty="dashed",lwd=3)
legend(x="topleft", legend=c("Literature 1", "Literature 2"),
        lty=1:2, cex=0.8)
dev.off()


png(file = "example_two_literatures_n.png",   # The directory you want to save the file in
    width = 400*2, # The width of the plot in inches
    height = 400)
par(mar=c(5,6,1,1)+.1)
plot(ns, pwrs1,type="line",cex.lab=1.8,lwd=3, ylab=TeX(r'(Average Power)'), xlab=TeX(r'(Sample Size of Every Experiment)'))
lines(ns,pwrs2,lty="dashed",lwd=3)
legend(x="topleft", legend=c("Literature 1", "Literature 2"),
       lty=1:2, cex=1.2)
dev.off()

#Maximum value of Delta
library(latex2exp)
pdf(file = "max_Delta.pdf",   # The directory you want to save the file in
    width = 4*1.6, # The width of the plot in inches
    height = 4)
cs <- 100:175 / 100
pwrs <- pnorm( 1.96*cs-1.96  )+pnorm(- 1.96*cs-1.96  )- pnorm( 1.96*1-1.96  )-pnorm(- 1.96*1-1.96  )
plot(cs^2,pwrs,cex.lab=1.5,lwd=3,type="line",ylab=TeX(r'(Max $\Delta_c$)'), xlab = TeX(r'($c^2$ (sample size increase factor) )') )
dev.off()

#VISUALIZE stretch and compress

ts <- -30:70/10

dev.off()





ts <- -30:60/10
lb1 <- 3
fun_pi0 <- function(t){
  #return( (abs(t) <= 0.25)+(t<=lb1+.5 & t>lb1) )
  return( (abs(t) <= 0.25)+(abs(t-3) <= 0.25) )
}
integrate(fun_pi0,-10,10)

lb2 <- lb1
#fun_pi1 <- function(t){
#  return( 2*(abs(t) <= 0.25/2)+2*(t<=lb2+.5*3/4 & t>lb2+.5/4) )
#}
fun_pi1 <- function(t){
  return( (abs(t) <= 0.4)-(abs(t)<=0.15)+(abs(t-3) <= 0.4)-(abs(t-3)<=0.15) )
}

c<- sqrt(2)
fun_K05_pi0 <- function(t){
  output <- t*0
  for (i in 1:length(t)){
    integrand <- function(tt){return(fun_pi0(tt)*dnorm((t[i]-tt)*c)*c)}
    output[i] <- integrate(integrand,-1,lb1+1)$value
  }
  return(output)
}
fun_K1_pi0 <- function(t){
  output <- t*0
  for (i in 1:length(t)){
    integrand <- function(tt){return(fun_pi0(tt)*dnorm((t[i]-tt)/1)/1)}
    output[i] <- integrate(integrand,-1,lb1+1)$value
  }
  return(output)
}
fun_K1_pi1 <- function(t){
  output <- t*0
  for (i in 1:length(t)){
    integrand <- function(tt){return(fun_pi1(tt)*dnorm((t[i]-tt)/1)/1)}
    output[i] <- integrate(integrand,-1,lb2+1)$value
  }
  return(output)
}
fun_Tc_pi1 <- function(t){
  output <- t*0
  for (i in 1:length(t)){
    integrand <- function(tt){return(fun_pi1(tt)*dnorm((t[i]-c*tt)/1)/1)}
    output[i] <- integrate(integrand,-1,lb2+1)$value
  }
  return(output)
}
fun_Tc_pi0 <- function(t){
  output <- t*0
  for (i in 1:length(t)){
    integrand <- function(tt){return(fun_pi0(tt)*dnorm((t[i]-c*tt)/1)/1)}
    output[i] <- integrate(integrand,-1,lb2+1)$value
  }
  return(output)
}
theta0 <- 0.5
fun_K1_pb_no_wt <- function(t){
  return(fun_K1_pi1(t)*((abs(t)<1.96)*theta0+(abs(t) >=1.96) ) )
}

fun_K1_pi0_pb <- function(t)
{
  return(fun_K1_pb_no_wt(t)/integrate(fun_K1_pb_no_wt,-6,6)$value )
}

ylims <- c(0,1.01)
xlims <- c(-2,6)
ca <- 2

lw <- 14
cl <- 4
cm <- 7
pdf(file = "conv_id_nine.pdf",   # The directory you want to save the file in
    width = 20*1.5, # The width of the plot in inches
    height = 20)
par(mfrow = c(2, 3),mar=c(8,8,8,8),oma=c(2,2,2,2),mgp = c(5, 1, 0))
plot(ts,fun_pi0(ts),type="l",ylim=ylims,xlim=xlims,cex.axis=ca,lwd=lw,cex.lab=cl,cex.main=cm,main=TeX(r'($\pi_1$)'),xlab=TeX(r'($h$)'),ylab="Density")
plot(ts,fun_K1_pi0(ts),ylim=ylims,xlim=xlims,type="l",cex.axis=ca,lwd=lw,cex.lab=cl,cex.main=cm,ylab="Density",main=TeX(r'($f_T$)'),xlab=TeX(r'($T$)'))
plot(ts,fun_Tc_pi0(ts),ylim=ylims,xlim=xlims,type="l",cex.axis=ca,lwd=lw,cex.lab=cl,cex.main=cm,ylab="Density",main=TeX(r'($f_{T_c}$)'),xlab=TeX(r'($T_c$)'))
plot(ts,fun_pi1(ts),type="l",ylim=ylims,xlim=xlims,cex.axis=ca,lwd=lw,cex.lab=cl,cex.main=cm,main=TeX(r'($\pi_2$)'),xlab=TeX(r'($h$)'),ylab="Density")
plot(ts,fun_K1_pi1(ts),ylim=ylims,xlim=xlims,type="l",cex.axis=ca,lwd=lw,cex.lab=cl,cex.main=cm,ylab="Density",main=TeX(r'($f_T$)'),xlab=TeX(r'($T$)'))
plot(ts,fun_Tc_pi1(ts),ylim=ylims,xlim=xlims,type="l",cex.axis=ca,lwd=lw,cex.lab=cl,cex.main=cm,ylab="Density",main=TeX(r'($f_{T_c}$)'),xlab=TeX(r'($T_c$)'))
dev.off()


ylims <- c(0,1.01)
pdf(file = "conv_id_four.pdf",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 20)
par(mfrow = c(2, 2),mar=c(8,8,8,8),oma=c(2,2,2,2),mgp = c(5, 1, 0))
plot(ts,fun_pi0(ts),type="l",ylim=ylims,xlim=xlims,lwd=lw,cex.axis=ca,cex.lab=cl,cex.main=cm,main=TeX(r'($\pi_1$)'),xlab=TeX(r'($h$)'),ylab="Density")
plot(ts,fun_K1_pi0(ts),ylim=ylims,xlim=xlims,type="l",cex.axis=ca,lwd=lw,cex.lab=cl,cex.main=cm,ylab="Density",main=TeX(r'($f_T$)'),xlab=TeX(r'($T$)'))
#plot(ts,fun_Tc_pi0(ts),ylim=ylims,xlim=xlims,type="l",lwd=lw,cex.lab=cl,cex.main=cm,ylab="Density",main=TeX(r'($f_{T_c}$)'),xlab=TeX(r'($T_c=ch_1+Z$)'))
plot(ts,fun_pi1(ts),type="l",ylim=ylims,xlim=xlims,lwd=lw,cex.axis=ca,cex.lab=cl,cex.main=cm,main=TeX(r'($\pi_2$)'),xlab=TeX(r'($h$)'),ylab="Density")
plot(ts,fun_K1_pi1(ts),ylim=ylims,xlim=xlims,type="l",cex.axis=ca,lwd=lw,cex.lab=cl,cex.main=cm,ylab="Density",main=TeX(r'($f_T$)'),xlab=TeX(r'($T$)'))
#plot(ts,fun_Tc_pi1(ts),ylim=ylims,xlim=xlims,type="l",lwd=lw,cex.lab=cl,cex.main=cm,ylab="Density",main=TeX(r'($f_{T_c}$)'),xlab=TeX(r'($T_c=ch_1+Z$)'))
dev.off()

pdf(file = "conv_id_three.pdf",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 8)
cm <- 5
par(mfrow = c(1, 3),mar=c(8,8,8,8),oma=c(3,3,3,3),mgp = c(5, 1, 0))
plot(ts,fun_pi0(ts),type="l",ylim=ylims,xlim=xlims,lwd=lw,cex.axis=ca,cex.lab=cl,cex.main=cm,main=TeX(r'($\pi_0$)'),xlab=TeX(r'($h$)'),ylab="Density")
plot(ts,fun_K1_pi0(ts),ylim=ylims,xlim=xlims,type="l",lwd=lw,cex.axis=ca,cex.lab=cl,cex.main=cm,ylab="Density",main=TeX(r'($f_T$ (PDF of $T))'),xlab=TeX(r'($T$)'))
plot(ts,fun_Tc_pi0(ts),ylim=ylims,xlim=xlims,type="l",lwd=lw,cex.axis=ca,cex.lab=cl,cex.main=cm,ylab="Density",main=TeX(r'($f_{T_c}$ (PDF of $T_c$))'),xlab=TeX(r'($T_c$)'))
#plot(ts,fun_pi1(ts),type="l",ylim=ylims,xlim=xlims,lwd=lw,cex.lab=cl,cex.main=cm,main=TeX(r'($\pi_2$)'),xlab=TeX(r'($h_2$)'),ylab="Density")
#plot(ts,fun_K1_pi1(ts),ylim=ylims,xlim=xlims,type="l",lwd=lw,cex.lab=cl,cex.main=cm,ylab="Density",main=TeX(r'($f_T$)'),xlab=TeX(r'($T=h_2+Z$)'))
#plot(ts,fun_Tc_pi1(ts),ylim=ylims,xlim=xlims,type="l",lwd=lw,cex.lab=cl,cex.main=cm,ylab="Density",main=TeX(r'($f_{T_c}$)'),xlab=TeX(r'($T_c=ch_1+Z$)'))
dev.off()


pdf(file = "conv_id_two_alt.pdf",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 8)
cm <- 5
par(mfrow = c(1, 3),mar=c(8,8,8,8),oma=c(3,3,3,3),mgp = c(5, 1, 0))
plot(ts,fun_pi0(ts),type="l",ylim=ylims,xlim=xlims,lwd=lw,cex.axis=ca,cex.lab=cl,cex.main=cm,main=TeX(r'($\pi_0$)'),xlab=TeX(r'($h$)'),ylab="Density")
plot(ts,fun_K1_pi0(ts),ylim=ylims,xlim=xlims,type="l",lwd=lw,cex.axis=ca,cex.lab=cl,cex.main=cm,ylab="Density",main=TeX(r'($f_T$ (PDF of $T))'),xlab=TeX(r'($T$)'))
dev.off()



pdf(file = "pb_three.pdf",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 8)
cm <- 5
ylims <- c(0,1)
par(mfrow = c(1, 3),mar=c(8,8,8,8),oma=c(3,3,3,3),mgp = c(5, 1, 0))
plot(ts,fun_pi0(ts),type="l",ylim=ylims,xlim=xlims,lwd=lw,cex.axis=ca,cex.lab=cl,cex.main=cm,main=TeX(r'($\pi_0$)'),xlab=TeX(r'($h$)'),ylab="Density")
plot(ts,fun_K1_pi0(ts),ylim=ylims,xlim=xlims,type="l",lwd=lw,cex.axis=ca,cex.lab=cl,cex.main=cm,ylab="Density",main=TeX(r'($f_T$ (no pb))'),xlab=TeX(r'($T$)'))
plot(ts,fun_K1_pi0_pb(ts),ylim=ylims,xlim=xlims,type="l",lwd=lw,cex.axis=ca,cex.lab=cl,cex.main=cm,ylab="Density",main=TeX(r'($f_{T}$ (with pb))'),xlab=TeX(r'($T$)'))
#plot(ts,fun_pi1(ts),type="l",ylim=ylims,xlim=xlims,lwd=lw,cex.lab=cl,cex.main=cm,main=TeX(r'($\pi_2$)'),xlab=TeX(r'($h_2$)'),ylab="Density")
#plot(ts,fun_K1_pi1(ts),ylim=ylims,xlim=xlims,type="l",lwd=lw,cex.lab=cl,cex.main=cm,ylab="Density",main=TeX(r'($f_T$)'),xlab=TeX(r'($T=h_2+Z$)'))
#plot(ts,fun_Tc_pi1(ts),ylim=ylims,xlim=xlims,type="l",lwd=lw,cex.lab=cl,cex.main=cm,ylab="Density",main=TeX(r'($f_{T_c}$)'),xlab=TeX(r'($T_c=ch_1+Z$)'))
dev.off()



pdf(file = "conv_id_two.pdf",   # The directory you want to save the file in
    width = 14, # The width of the plot in inches
    height = 8)
cm <- 2
ca <- 2
cl <- 2
par(mfrow = c(1, 2),mar=c(8,8,8,8),oma=c(3,3,3,3),mgp = c(5, 1, 0))
plot(ts,fun_pi0(ts),type="l",ylim=ylims,xlim=xlims,lwd=lw,cex.axis=ca,cex.lab=cl,cex.main=cm,main=TeX(r'($\pi_0$)'),xlab=TeX(r'($h$)'),ylab="Density")
plot(ts,fun_K1_pi0(ts),ylim=ylims,xlim=xlims,type="l",lwd=lw,cex.axis=ca,cex.lab=cl,cex.main=cm,ylab="Density",main=TeX(r'($f_T$ (PDF of $T))'),xlab=TeX(r'($T$)'))
dev.off()



#Check the above graphs with sims
sims <- 1000000
plot(density(( rnorm(2*sims)+c(c*runif(sims,min=-0.25,max=0.25),c*runif(sims,min=lb1,max=lb1+.5) ) )),ylim=c(0,1),xlim=c(-2,6))
plot(density(( rnorm(2*sims)/c+c(runif(sims,min=-0.25,max=0.25),runif(sims,min=lb1,max=lb1+.5) ) )),ylim=c(0,1),xlim=c(-3,5))
plot(density(( rnorm(2*sims)*0+c(runif(sims,min=-0.25,max=0.25),runif(sims,min=lb1,max=lb1+.5) ) )),ylim=c(0,1),xlim=c(-3,5))

source("C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/code/functions/Functions_TFUL.R")

pdf(file = "hermite.pdf",   # The directory you want to save the file in
    width = 4*1.6*1.5, # The width of the plot in inches
    height = 4)
ts <- -300:300 / 100
ylims <- c(-4,4)
lw <- 5
cl <- 1.2
plot(ts, hermite_general(ts,1,1),cex.lab=cl,main=TeX(r'(Hermite Polynomials: $\chi_j(t)$)'),ylim=ylims,type="l",lwd=lw,xlab="t",ylab=TeX(r'($\chi_j(t)$)'))
lines(ts, hermite_general(ts,3,1),type="l",lwd=lw,lty=2,col="blue")
#lines(ts, hermite_general(ts,3,1),type="l",lwd=lw,lty=3,col="green")
lines(ts, hermite_general(ts,10,1),type="l",lwd=lw,lty=4,col="red")
legend(x="bottomright", legend=c("j=1", "j=3","j=10"),
       col=c("black","blue","red"),lty=c(1,2,4), cex=0.8)
dev.off()


