# Testing for Underpowered Literatures
# Stefan Faridani
# February 26, 2024
# See README for replication instructions
# This script is for transparency purposes only. 
# This code will likely change prior to publication of the article.
# Run Main.R instead of this script.

#Replace root with appropriate path
setwd(root)

#load in functions
source("code/functions/Functions_TFUL.R")



#Many Labs

setwd(paste0(root,"/data"))
filename <- "ML-_Summary_Statistics.xlsx"
data_in = readxl::read_excel(filename,sheet="Anchoring1")
cpower_Anchoring1<- compute_cpower(data_in$Site, data_in$`Mean (High)`, data_in$`Mean (Low)`, data_in$`SD (High)`, data_in$`SD (Low)`, data_in$`N (High)`, data_in$`N (Low)`)
cpower_Anchoring1
powers <- cpower_Anchoring1
ML_tscores <- data_in$`t (unequal var)`
ML_sites <- data_in$Site
ML_treatments <- rep("Anchoring1",length(data_in$Site))

data_in = readxl::read_excel(filename,sheet="Anchoring2")
cpower_Anchoring2 <- compute_cpower(data_in$Site, data_in$`Mean (High)`, data_in$`Mean (Low)`, data_in$`SD (High)`, data_in$`SD (Low)`, data_in$`N (High)`, data_in$`N (Low)`)
cpower_Anchoring2
powers <- rbind(powers,cpower_Anchoring2)
ML_tscores <- c(ML_tscores,data_in$`t (unequal var)`)
ML_sites <- c(ML_sites,data_in$Site)
ML_treatments <- c(ML_treatments,rep("Anchoring2",length(data_in$Site)))

data_in = readxl::read_excel(filename,sheet="Anchoring3")
cpower_Anchoring3 <- compute_cpower(data_in$Site, data_in$`Mean (High)`, data_in$`Mean (Low)`, data_in$`SD (High)`, data_in$`SD (Low)`, data_in$`N (High)`, data_in$`N (Low)`)
cpower_Anchoring3
powers <- rbind(powers,cpower_Anchoring3)
ML_tscores <- c(ML_tscores,data_in$`t (unequal var)`)
ML_sites <- c(ML_sites,data_in$Site)
ML_treatments <- c(ML_treatments,rep("Anchoring3",length(data_in$Site)))

data_in = readxl::read_excel(filename,sheet="Anchoring4")
cpower_Anchoring4 <- compute_cpower(data_in$Site, data_in$`Mean (High)`, data_in$`Mean (Low)`, data_in$`SD (High)`, data_in$`SD (Low)`, data_in$`N (High)`, data_in$`N (Low)`)
cpower_Anchoring4
powers <- rbind(powers,cpower_Anchoring4)
ML_tscores <- c(ML_tscores,data_in$`t (unequal var)`)
ML_sites <- c(ML_sites,data_in$Site)
ML_treatments <- c(ML_treatments,rep("Anchoring4",length(data_in$Site)))

data_in = readxl::read_excel(filename,sheet="Gambler's Fallacy")
cpower_Gambler <- compute_cpower(data_in$Site, data_in$`Mean (Three6)`, data_in$`Mean (Two6)`, data_in$`SD (Three6)`, data_in$`SD (Two6)`, data_in$`N (Three6)`, data_in$`N (Two6)`)
cpower_Gambler
powers <- rbind(powers,cpower_Gambler)
ML_tscores <- c(ML_tscores,data_in$`t (unequal var)`)
ML_sites <- c(ML_sites,data_in$Site)
ML_treatments <- c(ML_treatments,rep("Gambler",length(data_in$Site)))

data_in = readxl::read_excel(filename,sheet="Money Priming")
cpower_Money <- compute_cpower(data_in$Site, data_in$`Mean (Money Primg)`, data_in$`Mean (Control)`, data_in$`SD (Money Primg)`, data_in$`SD (Control)`, data_in$`N (Money Priming)`, data_in$`N (Control)`)
cpower_Money
powers <- rbind(powers,cpower_Money)
ML_tscores <- c(ML_tscores,data_in$`t (unequal var)`)
ML_sites <- c(ML_sites,data_in$Site)
ML_treatments <- c(ML_treatments,rep("Money",length(data_in$Site)))

data_in = readxl::read_excel(filename,sheet="Imagined Contact")
cpower_Imagined <- compute_cpower(data_in$Site, data_in$`Mean (Contact)`, data_in$`Mean (Control)`, data_in$`SD (Contact)`, data_in$`SD (Control)`, data_in$`N (Contact)`, data_in$`N (Control)`)
cpower_Imagined
powers <- rbind(powers,cpower_Imagined)
ML_tscores <- c(ML_tscores,data_in$`t (unequal var)`)
ML_sites <- c(ML_sites,data_in$Site)
ML_treatments <- c(ML_treatments,rep("Imagined",length(data_in$Site)))

data_in = readxl::read_excel(filename,sheet="Sunk Costs")
cpower_Sunk <- compute_cpower(data_in$Site, data_in$`Mean (Paid)`, data_in$`Mean (Free)`, data_in$`SD (Paid)`, data_in$`SD (Free)`, data_in$`N (Paid)`, data_in$`N (Free)`)
cpower_Sunk
powers <- rbind(powers,cpower_Sunk)
ML_tscores <- c(ML_tscores,data_in$`t (unequal var)`)
ML_sites <- c(ML_sites,data_in$Site)
ML_treatments <- c(ML_treatments,rep("Sunk",length(data_in$Site)))

data_in = readxl::read_excel(filename,sheet="Quote Attribution")
cpower_Quote <- compute_cpower(data_in$Site, data_in$`Mean (Liked)`, data_in$`Mean (Disliked)`, data_in$`SD (Liked)`, data_in$`SD (Disliked)`, data_in$`N (Liked)`, data_in$`N (Disliked)`)
cpower_Quote
powers <- rbind(powers,cpower_Quote)
ML_tscores <- c(ML_tscores,data_in$`t (unequal var)`)
ML_sites <- c(ML_sites,data_in$Site)
ML_treatments <- c(ML_treatments,rep("Quote",length(data_in$Site)))

data_in = readxl::read_excel(filename,sheet="Flag Priming")
cpower_Flag <- compute_cpower(data_in$Site, data_in$`Mean (Flag Prime)`, data_in$`Mean (Control)`, data_in$`SD (Flag Prime)`, data_in$`SD (Control)`, data_in$`N (Flag Prime)`, data_in$`N (Control)`)
cpower_Flag
powers <- rbind(powers,cpower_Flag)
ML_tscores <- c(ML_tscores,data_in$`t (unequal var)`)
ML_sites <- c(ML_sites,data_in$Site)
ML_treatments <- c(ML_treatments,rep("Flag",length(data_in$Site)))

data_in = readxl::read_excel(filename,sheet="Math_Art Gender")
data_in <- data_in[c(1:22,24:28),]
data_in$`t (unequal var)` <- as.numeric(data_in$`t (unequal var)`)
cpower_Math <- compute_cpower(data_in$Site, data_in$`Mean (Female)`, data_in$`Mean (Male)`, data_in$`SD (Female)`, data_in$`SD (Male)`, data_in$`N (Female)`, data_in$`N (Male)`)
cpower_Math
powers <- rbind(powers,cpower_Math)
ML_tscores <- c(ML_tscores,data_in$`t (unequal var)`)
ML_sites <- c(ML_sites,data_in$Site)
ML_treatments <- c(ML_treatments,rep("Math",length(data_in$Site)))


ind <- ML_sites != "Overall for US participants:" & ML_sites !="Overall:" & ML_sites !="Mean across samples:" & ML_sites !="Overall (sum of samples)" 
C<- 2
D <- 1e-4
sigma_Y<-1
cv <- 1.96 #2.575 #1.96
c <- sqrt(2)

J_ML_sites     <- log(D*(length(unique(ML_sites[ind]) ))^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_ML_sites   <- C*(length(unique(ML_sites[ind]) ))^(-1/3)
ML_by_sites <- estimator(ML_tscores[ind],J=J_ML_sites,cv=cv,c=c,sigma_Y=1,bandwidth=eps_ML_sites,studies = ML_sites[ind],studies2= ML_treatments[ind])

J_ML_tscores     <- log(D*(length(ML_tscores[ind] ))^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_ML_tscores   <- C*(length(ML_tscores[ind]  ))^(-1/3)
ML_by_tscores <- estimator(ML_tscores[ind],J=J_ML_tscores,cv=cv,c=c,sigma_Y=1,bandwidth=eps_ML_tscores,studies = ML_sites[ind],studies2= ML_treatments[ind])

#robustness check
D <- 3e-4
C <- 1
J_ML_sites_r     <- log(D*(length(unique(ML_sites[ind]) ))^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_ML_sites_r   <- C*(length(unique(ML_sites[ind]) ))^(-1/3)
ML_by_sites_r <- estimator(ML_tscores[ind],J=J_ML_sites_r,cv=cv,c=c,sigma_Y=1,bandwidth=eps_ML_sites_r,studies = ML_sites[ind],studies2= ML_treatments[ind])

J_ML_tscores_r     <- log(D*(length(ML_tscores[ind] ))^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_ML_tscores_r   <- C*(length(ML_tscores[ind]  ))^(-1/3)
ML_by_tscores_r <- estimator(ML_tscores[ind],J=J_ML_tscores_r,cv=cv,c=c,sigma_Y=1,bandwidth=eps_ML_tscores_r,studies = ML_sites[ind],studies2= ML_treatments[ind])

mean(abs(ML_tscores[ind])>1.96)

#Alternative estimate assuming that true treatment effects are identical within experiments and across labs
deltas_by_intervention <- powers[,3]-powers[,2]
delta_hat_alt <- mean(deltas_by_intervention)
ses <- powers[,dim(powers)[2]]

#variance of alternative estimator assuming experiments within sites are independent
var <- 0
for (i in 1:length(ses)){
  for (j in 1:length(ses)){
    var <- var + ses[i]*ses[j]*(j==i)
  }
}
var <- var / length(ses)

#variance of alternative estimator assuming experiments within sites are perfectly dependent
var_cluster <- 0
for (i in 1:length(ses)){
  for (j in 1:length(ses)){
    var_cluster <- var_cluster + ses[i]*ses[j]
  }
}
var_cluster <- var_cluster / length(ses)

#Display results
print(paste0("Delta_hat_ALT: ", delta_hat_alt,", SE: ", sqrt(var),", SE Cluster:",sqrt(var_cluster) ))
table_ML <- cbind(ML_by_tscores,ML_by_sites,ML_by_tscores_r,ML_by_sites_r)
print(table_ML)
setwd(paste0(root,"/output/tables"))
write.csv(table_ML, file = 'results_ML.csv')
write.table(paste0("Delta_hat_ALT: ", mean(deltas),", SE: ", sqrt(var),", SE Cluster:",sqrt(var_cluster) ),file='results_ML_alt.txt')

#Compute ML power gain over different cs
cs <- sqrt(1+(0:10)/10)
deltas <- 0*cs
sd_deltas <- 0*cs
for (cc in 1:length(cs))
{
  print(paste0("It ", cc, " of ", length(cs)))
  
  ML<- estimator(ML_tscores[ind],J=J_ML_tscores,cv=cv,c=cs[cc],sigma_Y=1,bandwidth=eps_ML_tscores,studies = ML_sites[ind],studies2= ML_treatments[ind])
  deltas[cc] <- ML$deltahat
  sd_deltas[cc] <- ML$sd_delta

  print(cbind(cs[1:cc]^2, deltas[1:cc], sd_deltas[1:cc]))
}
setwd(paste0(root))
makeciplot_single_gg(cs,deltas,sd_deltas,"Many Labs Replications")
ggsave("output/figures/Single_Plot_ML.pdf")

makeciplot_double_gg(cs,deltas,deltas_rcts,sd_deltas,sd_deltas_rcts,"Many Labs","RCTs","")
ggsave("output/figures/ML_vs_RCTs.pdf")