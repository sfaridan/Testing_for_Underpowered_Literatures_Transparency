rm(list = ls())

#load in functions
source("C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/code/functions/Functions_TFUL.R")

#Load in the data
setwd("C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/data")
MM_DID = read.csv('MM_DID.csv')
MM = read.csv('MM_RCTs.csv')
MM_DD = read.csv('MM_DD.csv')
MM_IV = read.csv('MM_IV.csv')
MM_other = rbind(MM_DD,MM_IV,MM_DID)
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


table_values(out_RCT_articles)
table_values(out_RCT_tscores)
table_values(out_other_articles)
table_values(out_other_tscores)
c(J_RCT_articles, J_other_articles, J_RCT_tscores, J_other_tscores)
c(eps_RCT_articles, eps_other_articles, eps_RCT_articles, eps_other_tscores)

#test for equality
2*(1-pnorm((out_RCT_articles$deltahat-out_other_articles$deltahat)/sqrt(out_RCT_articles$varest_delta+out_other_articles$varest_delta  )))
2*(1-pnorm((out_RCT_tscores$deltahat-out_other_tscores$deltahat)/sqrt(out_RCT_tscores$varest_delta+out_other_tscores$varest_delta  )))

#going from 5% to 1% for c=sqrt(2): slight increase in sensitivity (9.6%) but still sig more for non-rcts
# c=sqrt(1.5): same story. Sensitivity 3.6% and still significantly more for non-rcts
#comparisons with Many Labs hold up for both
