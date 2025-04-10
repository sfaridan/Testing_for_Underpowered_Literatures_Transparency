# Testing for Underpowered Literatures
# Stefan Faridani
# February 26, 2024
# See README for replication instructions
# This script is for transparency purposes only. 
# This code will likely change prior to publication of the article.
# Run this on its own (to avoid interference from GGplot)


#generate t-scores
setwd(paste0(root,"/output/figures"))

n <- 50000
theta <- 0.5
xlims <- c(-5,5)
ylims <- c(0,0.4)
cexlab <- 1.4
hs <- c(-1+rnorm(n/2),2+0.6*rnorm(n/2))
hs_big <- c(-1+rnorm(100000),2+0.6*rnorm(100000))
Ts <- hs + rnorm(n)
include <- abs(Ts)>1.96 | runif(n)<theta
Ts_pb <- Ts[ include ]

breakpoints <- ((-9*3):(10*3))/3 #(-90:100)/10


pdf(file = "Slides_Tscores_pb.pdf",   # The directory you want to save the file in
)
hist(Ts_pb,xlim =xlims,ylim=ylims,
     breaks=breakpoints,xlab="t-scores",
     main="Sample of T-scores",
     probability=TRUE,cex.lab=cexlab)
dev.off()

pdf(file = "Slides_Tscores_nopb.pdf",   # The directory you want to save the file in
)
hist(Ts,xlim =xlims,ylim=ylims,
     breaks=breakpoints,xlab="t-scores",
     main="Pub Bias Removed",probability=TRUE,cex.lab=cexlab)
dev.off()

pdf(file = "Slides_fit.pdf",   # The directory you want to save the file in
)
hist(Ts,xlim =xlims,ylim=ylims,
     breaks=breakpoints,xlab="t-scores",
     main="",probability=TRUE,cex.lab=cexlab)
lines(density(hs_big), # density plot
      lwd = 3, # thickness of line
      col = "blue")
legend(x="topright",y="right",legend=c("True Effects"),
       col=c("blue"), lty=1, cex=1.3)
dev.off()

pdf(file = "Slides_True_Effects_Only.pdf",   # The directory you want to save the file in
    )
plot(density(hs_big),xlim =xlims,ylim=ylims,
     lwd = 5, # thickness of line
     col = "blue",xlab="",main="",cex.lab=cexlab)
legend(x="topright",y="right",legend=c("True Effects"),
       col=c("blue"), lty=1, cex=1.3)
dev.off()


pdf(file = "Slides_True_Effects_T_counter.pdf",   # The directory you want to save the file in
)
plot(density(hs_big),xlim =xlims,ylim=ylims,
     lwd = 5, # thickness of line
     col = "blue",xlab="",main="",cex.lab=cexlab)
lines(density(hs_big*sqrt(2)+rnorm(length(hs_big))), # density plot
      lwd = 5, # thickness of line
      col = "red",xlab="",main="")
legend(x="topright",y="right",legend=c("True Effects", "T-scores (double n)"),
       col=c("blue", "red"), lty=1, cex=1.3)
dev.off()

pdf(file = "Slides_all.pdf",   # The directory you want to save the file in
)
hist(Ts,xlim =xlims,ylim=ylims,
     breaks=breakpoints,xlab="t-scores",
     main="",probability=TRUE,cex.lab=cexlab)
lines(density(hs_big), # density plot
      lwd = 3, # thickness of line
      col = "blue")
lines(density(hs_big*sqrt(2)+rnorm(length(hs_big))), # density plot
      lwd = 5, # thickness of line
      col = "red",xlab="",main="")
legend(x="topright",y="right",legend=c("True Effects", "T-scores (double n)"),
       col=c("blue", "red"),lty=c(1,1), cex=1.3)
dev.off()
