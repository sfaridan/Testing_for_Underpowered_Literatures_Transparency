ns <- (10:30) *10
pwrs1 <- pnorm(sqrt(ns)*.17 -1.65)
pwrs2 <- 0.5*pnorm(sqrt(ns)*.04 -1.65)+0.5*pnorm(sqrt(ns)*2 -1.65)

setwd(paste0(here(),"/R/Underpowered Literatures/output/figures" ))

 
plot(ns, pwrs1,type="line",cex.lab=1, ylab="Power", xlab="Sample size (n)")
lines(ns,pwrs2,lty="dashed")
legend(x="topleft", legend=c("Literature 1", "Literature 2"),
        lty=1:2, cex=0.8)

