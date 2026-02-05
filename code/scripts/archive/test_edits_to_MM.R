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
#root <- "C:/Users/sfaridani6/Documents/GitHub/Testing_for_Underpowered_Literatures_Transparency"
setwd(root)

#load in functions
source(paste0(root,"/code/functions/Functions_TFUL_new_inference.R"))
library(ggplot2)

#Load in the data
setwd(data)
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
J_RCT_tscores     <- log(D*(length(MM$t))^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_RCT_tscores <- C*(length(MM$t))^(-1/3)
out_RCT_tscores <- estimator(MM$t,J=J_RCT_tscores,cv=cv,c=c,sigma_Y=1,bandwidth=eps_RCT_tscores,studies = MM$unique_id,include_pb=TRUE)
out_RCT_tscores

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
#write.csv(table_econ, file = 'results_MM.csv')

#display results nicely
ses<- sqrt(c(out_RCT_tscores$varest_delta,out_other_tscores$varest_delta,out_RCT_articles$varest_delta,out_other_articles$varest_delta ))
deltas <- c(out_RCT_tscores$deltahat,out_other_tscores$deltahat,out_RCT_articles$deltahat,out_other_articles$deltahat )
citop <- deltas+1.96*ses
cibot <- deltas-1.96*ses
results_old <- -cbind(deltas, -ses,citop,cibot)


####################
####################

#Set the seed (for consistent de-rounding)
set.seed(1) 

#Replace root with appropriate path
#root <- "C:/Users/sfaridani6/Documents/GitHub/Testing_for_Underpowered_Literatures_Transparency"
setwd(root)

#load in functions
source(paste0(root,"/code/functions/Functions_TFUL_Feb_2026.R"))
library(ggplot2)

#Load in the data
setwd(data)
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
J_RCT_tscores     <- log(D*(length(MM$t))^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_RCT_tscores <- C*(length(MM$t))^(-1/3)
out_RCT_tscores <- estimator(MM$t,J=J_RCT_tscores,cv=cv,c=c,sigma_Y=1,bandwidth=eps_RCT_tscores,studies = MM$unique_id,include_pb=TRUE)
out_RCT_tscores

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
#write.csv(table_econ, file = 'results_MM.csv')

#display results nicely
ses<- sqrt(c(out_RCT_tscores$varest_delta,out_other_tscores$varest_delta,out_RCT_articles$varest_delta,out_other_articles$varest_delta ))
deltas <- c(out_RCT_tscores$deltahat,out_other_tscores$deltahat,out_RCT_articles$deltahat,out_other_articles$deltahat )
citop <- deltas+1.96*ses
cibot <- deltas-1.96*ses
results_new <- cbind(deltas, ses,cibot,citop)
round(results,3)
round(results_new,3)