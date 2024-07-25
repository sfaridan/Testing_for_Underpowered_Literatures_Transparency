# Testing for Underpowered Literatures
# Stefan Faridani
# February 26, 2024
# See README for replication instructions
# This script is for transparency purposes only. 
# This code will likely change prior to publication of the article.
# Run Main.R instead of this script.

#Set the seed (for consistent de-rounding)
set.seed(1) 

#Replace root with appropriate path
#root <- "C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures"
setwd(root)

#load in functions
source(paste0(root,"/code/functions/Functions_TFUL.R"))
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

#partition data into rcts and non-rcts
is_rct <- MM_data$method == "RCT"
is_not_rct <- MM_data$method == "RDD" | MM_data$method == "DID" | MM_data$method == "IV"
MM <- MM_data[is_rct,]
MM_other <-MM_data[is_not_rct,]
MM_articles <- length(unique (MM$article))
Other_articles <- length(unique (MM_other$article))

#Preferred tuning constants
C <- 2
D <- 1e-4
sigma_Y <- 1
cv <-  1.96 #2.575
c <- sqrt(2) #sqrt(1.5)

#Brodeur RCTs, scale tuning parameters with number of articles
J_RCT_articles     <- log(D*(MM_articles)^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_RCT_articles <- C*(MM_articles)^(-1/3)
out_RCT_articles <- estimator(MM$t,J=J_RCT_articles,cv=cv,c=c,sigma_Y=1,bandwidth=eps_RCT_articles,studies = MM$article,include_pb=1)

#Brodeur RCTs, scale tuning parameters with number of t-scores
J_RCT_tscores     <- log(D*(length(MM$t))^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_RCT_tscores <- C*(length(MM$t))^(-1/3)
out_RCT_tscores <- estimator(MM$t,J=J_RCT_tscores,cv=cv,c=c,sigma_Y=1,bandwidth=eps_RCT_tscores,studies = MM$article,include_pb=1)

#Brodeur DID, IV, and RDD, scale tuning parameters with number of articles
J_other_articles     <- log(D*(Other_articles)^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_other_articles <- C*(Other_articles)^(-1/3)
out_other_articles <- estimator(MM_other$t,J=J_other_articles,cv=cv,c=c,sigma_Y=1,bandwidth=eps_other_articles,studies = MM_other$article,include_pb=1)

#Brodeur DID, IV, and RDD, scale tuning parameters with number of t-scores
J_other_tscores     <- log(D*(length(MM_other$t))^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_other_tscores <- C*(length(MM_other$t))^(-1/3)
out_other_tscores <- estimator(MM_other$t,J=J_other_tscores,cv=cv,c=c,sigma_Y=1,bandwidth=eps_other_tscores,studies = MM_other$article,include_pb=1)

 

#display main table results
table_econ <- cbind(out_RCT_tscores,out_other_tscores,out_RCT_articles,out_other_articles)

#Save tables
setwd(paste0(root,"/output/tables"))
write.csv(table_econ, file = 'results_MM.csv')


#Make plots
cs <- sqrt(1+(0:10)/10)
deltas_rcts <- 0*cs
deltas_nonrcts <- 0*cs
sd_deltas_rcts <- 0*cs
sd_deltas_nonrcts <- 0*cs
for (cc in 1:length(cs))
{
  print(paste0("It ", cc, " of ", length(cs)))
  
  rcts<- estimator(MM$t,J=J_RCT_tscores,cv=cv,c=cs[cc],sigma_Y=1,bandwidth=eps_RCT_tscores,studies = MM$article,include_pb=1)
  deltas_rcts[cc] <- rcts$deltahat
  sd_deltas_rcts[cc] <- rcts$sd_delta
  
  nonrcts <- estimator(MM_other$t,J=J_other_tscores,cv=cv,c=cs[cc],sigma_Y=1,bandwidth=eps_other_tscores,studies = MM_other$article,include_pb=1)
  deltas_nonrcts[cc] <- nonrcts$deltahat
  sd_deltas_nonrcts[cc] <- nonrcts$sd_delta
 
  print(cbind(cs[1:cc]^2, deltas_rcts[1:cc], sd_deltas_rcts[1:cc],deltas_nonrcts[1:cc],  sd_deltas_nonrcts[1:cc]))
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

