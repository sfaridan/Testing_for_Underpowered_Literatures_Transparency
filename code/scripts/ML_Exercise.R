#Testing for Underpowered Literatures
#Stefan Faridani
#February 26, 2024
#See README for replication instructions

rm(list = ls())

#Replace root with appropriate path
root <- "C:/Users/stefa/OneDrive/Documents/R/Underpowered_Literatures_Transparency"
setwd(root)

#load in functions
source("code/functions/Functions_TFUL.R")

#Computes counterfactual power for a set of experiments all studying the same intervention
compute_cpower_unweighted <- function(sites, mean1, mean2, sd1, sd2, n1, n2,c=sqrt(2)){
  ind <- sites != "Overall for US participants:" & sites !="Overall:" & sites !="Mean across samples:" & sites !="Overall (sum of samples)" 
  es <- mean(mean1[ind]-mean2[ind])
  hs <- es / sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind])
  pwr_status_quo <- 1-pnorm(1.96-hs)+pnorm(-hs-1.96)
  pwr_c <- 1-pnorm(1.96-c*hs)+pnorm(-c*hs-1.96)
  rejs <- abs((mean1[ind]-mean2[ind])/sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind]))>1.96
  #uncertainty <- sqrt(sd1[1]^2/n1[1]+sd2[1]^2/n2[1]) / mean(sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind]))
  #print((mean1[ind]-mean2[ind])/sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind]))
  
  #compute standard error
  se_es <- sqrt(sum(sd1[ind]^2/n1[ind]) / length(sd1[ind])^2 +sum(sd2[ind]^2/n2[ind]) / length(sd2[ind])^2   )
  se_hs <- se_es / sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind]) 
  se_pwrs <- se_hs * (dnorm(1.96-c*hs)+dnorm(-c*hs-1.96))  #Taylor approximation
  se_pwr <- mean(se_pwrs) /sqrt(length(se_pwrs))
  print(paste0("se / es: ", se_es/es, ", se_pwr: ",se_pwr, ", delta: ",mean(pwr_status_quo)- mean(pwr_c) ))
  
  return(c(mean(rejs),mean(pwr_status_quo), mean(pwr_c), se_pwr    ))
}

compute_cpower <- function(sites, mean1, mean2, sd1, sd2, n1, n2,c=sqrt(2)){
  ind <- sites != "Overall for US participants:" & sites !="Overall:" & sites !="Mean across samples:" & sites !="Overall (sum of samples)" 
  weights <- (n1[ind]+n2[ind])/sum(n1[ind]+n2[ind])
  es <- sum((mean1[ind]-mean2[ind])*weights)
  hs <- es / sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind])
  pwr_status_quo <- 1-pnorm(1.96-hs)+pnorm(-hs-1.96)
  pwr_c <- 1-pnorm(1.96-c*hs)+pnorm(-c*hs-1.96)
  rejs <- abs((mean1[ind]-mean2[ind])/sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind]))>1.96
  #uncertainty <- sqrt(sd1[1]^2/n1[1]+sd2[1]^2/n2[1]) / mean(sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind]))
  #print((mean1[ind]-mean2[ind])/sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind]))
  
  #compute standard error
  se_es <- sqrt(sum(sd1[ind]^2/n1[ind]*weights^2) +sum(sd2[ind]^2/n2[ind]*weights^2)   )
  se_hs <- se_es / sqrt(sd1[ind]^2/n1[ind]+sd2[ind]^2/n2[ind]) 
  se_pwrs <- se_hs * (dnorm(1.96-c*hs)+dnorm(-c*hs-1.96))  #Taylor approximation
  se_pwr <- mean(se_pwrs) /sqrt(length(se_pwrs))
  print(paste0("se / es: ", se_es/es, ", se_pwr: ",se_pwr, ", delta: ",mean(pwr_status_quo)- mean(pwr_c) ))
  
  return(c(mean(rejs),mean(pwr_status_quo), mean(pwr_c), se_pwr    ))
}

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


deltas <- powers[,3]-powers[,2]
ses <- powers[,dim(powers)[2]]
powers <- cbind(powers,deltas)
powers
mean(deltas)
var <- 0
for (i in 1:length(ses)){
  for (j in 1:length(ses)){
    var <- var + ses[i]*ses[j]*(j==i)
  }
}
var <- var / length(ses)
print(paste0("Delta_hat_ALT: ", mean(deltas),", SE: ", sqrt(var) ))
#sqrt(mean(abs(powers[,1]-powers[,2])^2))

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


table_values(ML_by_sites)
table_values(ML_by_sites_r)
table_values(ML_by_tscores)
table_values(ML_by_tscores_r)
print(c(J_ML_sites, J_ML_sites_r, J_ML_tscores, J_ML_tscores_r))
print(c(eps_ML_sites, eps_ML_sites_r, eps_ML_tscores, eps_ML_tscores_r))
mean(abs(ML_tscores[ind])>1.96)
