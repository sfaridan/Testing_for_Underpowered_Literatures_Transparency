# Testing for Underpowered Literatures
# Stefan Faridani
# sfaridani6@gatech.edu
# April 10, 2025
# See README for replication instructions
# This script is for transparency purposes only. 
# This code will likely change prior to publication of the article.
# Run Main.R instead of this script.

#Set the seed (for consistent de-rounding)
set.seed(1) 

#Replace root with appropriate path
root <- "C:/Users/sfaridani6/Documents/GitHub/Testing_for_Underpowered_Literatures_Transparency"
setwd(root)

#load in functions
source(paste0(root,"/code/functions/Functions_TFUL_new_inference.R"))
library(ggplot2)

#Load in the data
setwd(paste0(root,"/data"))
MM_data <- haven::read_dta('MM data.dta')

#drop observations where the point estimate or standard error are missing (can't de-round) 
to_drop <-  is.na(MM_data$mu) |  is.na(MM_data$sd) # | MM_data$sd == 0  
to_drop[is.na(to_drop)] <- 1
sum(to_drop)
#MM_data<- MM_data[!to_drop,]

#de-round
MM_data$mu <- de_round(MM_data$mu)
MM_data$sd <- de_round(MM_data$sd)
MM_data$t_der <- MM_data$mu/MM_data$sd
MM_data$t[!is.na(MM_data$t_der)] <- MM_data$t_der[!is.na(MM_data$t_der)] #replace t with derounded whenever it exists

#de-round everything else
MM_data$t[is.na(MM_data$t_der)] <- de_round(MM_data$t[is.na(MM_data$t_der)])

#uniquely identify the articles (do not use the article variable)
MM_data$unique_id <- as.factor(MM_data$title)
length(unique(MM_data$unique_id))

#partition data into rcts and non-rcts
is_rct <- MM_data$method == "RCT"
is_not_rct <- MM_data$method == "RDD" | MM_data$method == "DID" | MM_data$method == "IV"
MM <- MM_data[is_rct,]
MM_other <-MM_data[is_not_rct,]
MM_articles <- length(unique (MM$unique_id))
Other_articles <- length(unique (MM_other$unique_id))

#Preferred tuning constants
C <- 2
D <- 1e-4
sigma_Y <- 1
cv <-  1.96 #2.575
c <- sqrt(2) #sqrt(1.5)

#Brodeur RCTs, scale tuning parameters with number of t-scores
#On Faridani's desktop, this block takes 2.2 seconds to run
tic()
J_RCT_tscores     <- log(D*(length(MM$t))^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_RCT_tscores <- C*(length(MM$t))^(-1/3)
out_RCT_tscores <- estimator(MM$t,J=J_RCT_tscores,cv=cv,c=c,sigma_Y=1,bandwidth=eps_RCT_tscores,studies = MM$unique_id,include_pb=TRUE)
out_RCT_tscores
toc()

#Brodeur DID, IV, and RDD, scale tuning parameters with number of t-scores
J_other_tscores     <- log(D*(length(MM_other$t))^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_other_tscores <- C*(length(MM_other$t))^(-1/3)
out_other_tscores <- estimator(MM_other$t,J=J_other_tscores,cv=cv,c=c,sigma_Y=1,bandwidth=eps_other_tscores,studies = MM_other$unique_id,include_pb=TRUE)
out_other_tscores

#Brodeur DID, IV, and RDD, scale tuning parameters with number of articles
J_other_articles     <- log(D*(Other_articles)^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_other_articles <- C*(Other_articles)^(-1/3)
out_other_articles <- estimator(MM_other$t,J=J_other_articles,cv=cv,c=c,sigma_Y=1,bandwidth=eps_other_articles,studies = MM_other$unique_id,include_pb=TRUE)
out_other_articles

#Brodeur RCTs, scale tuning parameters with number of articles
J_RCT_articles     <- log(D*(MM_articles)^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_RCT_articles <- C*(MM_articles)^(-1/3)
out_RCT_articles <- estimator(MM$t,J=J_RCT_articles,cv=cv,c=c,sigma_Y=1,bandwidth=eps_RCT_articles,studies = MM$unique_id,include_pb=TRUE)
out_RCT_articles

#display main table results
table_econ <- cbind(out_RCT_tscores,out_other_tscores,out_RCT_articles,out_other_articles)

#Save tables
setwd(paste0(root,"/output/tables"))
write.csv(table_econ, file = 'results_MM.csv')

#display results nicely
ses<- sqrt(c(out_RCT_tscores$varest_delta,out_other_tscores$varest_delta,out_RCT_articles$varest_delta,out_other_articles$varest_delta ))
deltas <- c(out_RCT_tscores$deltahat,out_other_tscores$deltahat,out_RCT_articles$deltahat,out_other_articles$deltahat )
citop <- deltas+1.96*ses
cibot <- deltas-1.96*ses
results <- -cbind(deltas, -ses,citop,cibot)
round(results,3)

#Make plots
cs <- sqrt(1+(0:10)/10)
deltas_rcts <- 0*cs
deltas_nonrcts <- 0*cs
sd_deltas_rcts <- 0*cs
sd_deltas_nonrcts <- 0*cs
for (cc in 2:length(cs)) #we already know that for c=1, delta=0
{
  tic()
  print(paste0("It ", cc, " of ", length(cs)))
  
  rcts<- estimator(MM$t,J=J_RCT_tscores,cv=cv,c=cs[cc],sigma_Y=1,bandwidth=eps_RCT_tscores,studies = MM$unique_id,include_pb=TRUE)
  deltas_rcts[cc] <- rcts$deltahat
  sd_deltas_rcts[cc] <- sqrt(rcts$varest_delta)
  
  nonrcts <- estimator(MM_other$t,J=J_other_tscores,cv=cv,c=cs[cc],sigma_Y=1,bandwidth=eps_other_tscores,studies = MM_other$unique_id,include_pb=TRUE)
  deltas_nonrcts[cc] <- nonrcts$deltahat
  sd_deltas_nonrcts[cc] <- sqrt(nonrcts$varest_delta)
  
  print(cbind(cs[1:cc]^2, deltas_rcts[1:cc], sd_deltas_rcts[1:cc],deltas_nonrcts[1:cc],  sd_deltas_nonrcts[1:cc]))
  write.table(cbind(deltas_rcts,sd_deltas_rcts,deltas_nonrcts,sd_deltas_nonrcts,cc),file='MM_plot_as_table.txt',sep=" ")
  toc()
}
setwd(paste0(root))
makeciplot_double_gg(cs,deltas_rcts,deltas_nonrcts,sd_deltas_rcts,sd_deltas_nonrcts, "RCTs","Quasi-Experiments","")
ggsave("output/figures/Two_Plots_MM.pdf")

makeciplot_double_gg(cs,deltas_rcts,-100*deltas_nonrcts,sd_deltas_rcts,0*sd_deltas_nonrcts, "RCTs","Quasi-Experiments","")
ggsave("output/figures/One_Plots_MM.jpg")

#print significance of differences over cs
setwd(paste0(root,"/output/tables"))
pval_difference <- 2*(1-pnorm(abs(deltas_rcts-deltas_nonrcts)/sqrt(sd_deltas_rcts^2+sd_deltas_nonrcts^2)))
print(pval_difference)
write.table(cbind(cs,pval_difference),file='pvals_rct_vs_nonrct_overcs.txt',sep=" ")



#zcurve
set.seed(1)
power_function <- function(mu, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha / 2)
  power <- pnorm(-z_alpha - mu) + 1 - pnorm(z_alpha - mu)
  return(power)
}
library(zcurve)

#Zcurve non-rcts
meta_df <- MM_other
ts <- as.vector(meta_df$t)
tstr <- paste("z= ",as.character(meta_df$t))
zd <- zcurve_data(tstr,id=as.vector(meta_df$article))
zc_non <- zcurve(ts,bootstrap=FALSE)
delta_zcurve_nonrct  <- (power_function(sqrt(2)*zc_non$fit$mu)-power_function(zc_non$fit$mu))%*%zc_non$fit$weights
delta_zcurve_nonrct
#ses for zcurve for nonrcts
boots<-100
dboots <- rep(0,boots)
for (i in 1:boots) {
  #tosamp <- sample(1:(dim(MM)[1]),dim(MM)[1],replace = TRUE)
  cluster_ids <- unique(meta_df$unique_id)
  sampled_ids <- sample(cluster_ids, size = length(cluster_ids), replace = TRUE)
  resampled_list <- vector("list", length(sampled_ids))
  for (j in seq_along(sampled_ids)) {
    rows <- meta_df[meta_df$unique_id == sampled_ids[j], ]
    rows$replicate <- j  # Track which bootstrap replicate this is
    resampled_list[[j]] <- rows
  }
  resampled_df <- do.call(rbind, resampled_list)
  tstr <- paste("z= ",as.character(resampled_df$t))
  #zd <- zcurve_data(tstr,id=as.vector(resampled_df$article))
  zc_boot <- zcurve(as.vector(resampled_df$t) ,bootstrap=FALSE)
  #zc_boot <- zcurve_clustered(zd,bootstrap=100,method = "b")
  mus_boot <- zc_boot$fit$mu 
  dboots[i] <- (power_function(sqrt(2)*mus_boot)-power_function(mus_boot)) %*%zc_boot$fit$weights
  print(c(i,sd(dboots[1:i])))
}
se_nonrct <- sd(dboots)
c(delta_zcurve_nonrct,se_nonrct)

#Zcurve rcts
set.seed(1)
meta_df <- MM
ts <- as.vector(meta_df$t)
tstr <- paste("z= ",as.character(meta_df$t))
zd <- zcurve_data(tstr,id=as.vector(meta_df$article))
zc_rct <- zcurve(ts,bootstrap=FALSE)
delta_zcurve_rct  <- (power_function(sqrt(2)*zc_rct$fit$mu)-power_function(zc_rct$fit$mu))%*%zc_rct$fit$weights
delta_zcurve_rct
#ses for zcurve for nonrcts
boots<-100
dboots <- rep(0,boots)
for (i in 1:boots) {
  #tosamp <- sample(1:(dim(MM)[1]),dim(MM)[1],replace = TRUE)
  
  
  cluster_ids <- unique(meta_df$unique_id)
  sampled_ids <- sample(cluster_ids, size = length(cluster_ids), replace = TRUE)
  resampled_list <- vector("list", length(sampled_ids))
  for (j in seq_along(sampled_ids)) {
    rows <- meta_df[meta_df$unique_id == sampled_ids[j], ]
    rows$replicate <- j  # Track which bootstrap replicate this is
    resampled_list[[j]] <- rows
  }
  resampled_df <- do.call(rbind, resampled_list)
  tstr <- paste("z= ",as.character(resampled_df$t))
  zd <- zcurve_data(tstr,id=as.vector(resampled_df$article))
  zc_boot <- zcurve(as.vector(resampled_df$t) ,bootstrap=FALSE)
  #zc_boot <- zcurve_clustered(zd,bootstrap=100,method = "b")
  mus_boot <- zc_boot$fit$mu 
  dboots[i] <- (power_function(sqrt(2)*mus_boot)-power_function(mus_boot)) %*%zc_boot$fit$weights
  print(c(i,sd(dboots[1:i])))
}
se_rct <- sd(dboots)
c(delta_zcurve_rct,se_rct)

# full results
table_zcurve<-cbind(c(delta_zcurve_rct,se_rct),c(delta_zcurve_nonrct,se_nonrct) )
write.csv(table_zcurve, file = 'results_zcurve.csv')

#clustered results
#> c(delta,se_cluster)
#[1] 0.14107299 0.01518378
#> delta +1.96*c(-se_cluster,se_cluster)
#[1] 0.1113128 0.1708332
