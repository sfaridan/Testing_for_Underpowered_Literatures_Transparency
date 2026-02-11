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
setwd(root)

#load in functions
source(paste0(root,"/code/functions/Functions_TFUL.R"))

#Load in the data
setwd(data)
MM_data <- haven::read_dta('MM data.dta')

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
is_rct         <- MM_data$method == "RCT"
is_not_rct     <- MM_data$method == "RDD" | MM_data$method == "DID" | MM_data$method == "IV"
is_rct[is.na(is_rct)]         <- FALSE
is_not_rct[is.na(is_not_rct)] <- FALSE
MM             <- MM_data[is_rct,]
MM_other       <- MM_data[is_not_rct,]
MM_articles    <- length(unique (MM$unique_id))
Other_articles <- length(unique (MM_other$unique_id))

#Preferred tuning constants
C       <- 2
D       <- 1e-4
sigma_Y <- 1
cv      <-  1.96  
c       <- sqrt(2)  

#Brodeur RCTs, scale tuning parameters with number of t-scores
#On Faridani's desktop, this block takes 2.2 seconds to run
J_RCT_tscores   <- log(D*(length(MM$t))^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_RCT_tscores <- C*(length(MM$t))^(-1/3)
out_RCT_tscores <- estimator(MM$t,J=J_RCT_tscores,cv=cv,c=c,sigma_Y=1,bandwidth=eps_RCT_tscores,studies = MM$unique_id,include_pb=TRUE)
out_RCT_tscores # main estimate

#Brodeur DID, IV, and RDD, scale tuning parameters with number of t-scores
J_other_tscores     <- log(D*(length(MM_other$t))^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_other_tscores <- C*(length(MM_other$t))^(-1/3)
out_other_tscores <- estimator(MM_other$t,J=J_other_tscores,cv=cv,c=c,sigma_Y=1,bandwidth=eps_other_tscores,studies = MM_other$unique_id,include_pb=TRUE)
out_other_tscores

#Brodeur DID, IV, and RDD, scale tuning parameters with number of articles
J_other_articles   <- log(D*(Other_articles)^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_other_articles <- C*(Other_articles)^(-1/3)
out_other_articles <- estimator(MM_other$t,J=J_other_articles,cv=cv,c=c,sigma_Y=1,bandwidth=eps_other_articles,studies = MM_other$unique_id,include_pb=TRUE)
out_other_articles

#Brodeur RCTs, scale tuning parameters with number of articles
J_RCT_articles   <- log(D*(MM_articles)^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps_RCT_articles <- C*(MM_articles)^(-1/3)
out_RCT_articles <- estimator(MM$t,J=J_RCT_articles,cv=cv,c=c,sigma_Y=1,bandwidth=eps_RCT_articles,studies = MM$unique_id,include_pb=TRUE)
out_RCT_articles

#print results
ses<- sqrt(c(out_RCT_tscores$varest_delta,out_other_tscores$varest_delta,out_RCT_articles$varest_delta,out_other_articles$varest_delta ))
deltas <- c(out_RCT_tscores$deltahat,out_other_tscores$deltahat,out_RCT_articles$deltahat,out_other_articles$deltahat )
citop <- deltas+1.96*ses
cibot <- deltas-1.96*ses
results <- -cbind(deltas, ses,cibot,citop)
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
ggsave("output/figures/Two_Plots_MM.pdf",width=8,height=5)

#makeciplot_double_gg(cs,deltas_rcts,-100*deltas_nonrcts,sd_deltas_rcts,0*sd_deltas_nonrcts, "RCTs","Quasi-Experiments","")
#ggsave("output/figures/One_Plots_MM.jpg")

#print significance of differences over cs
setwd(paste0(root,"/output/tables"))
pval_difference <- 2*(1-pnorm(abs(deltas_rcts-deltas_nonrcts)/sqrt(sd_deltas_rcts^2+sd_deltas_nonrcts^2)))
print(pval_difference)
write.table(cbind(cs,pval_difference),file='pvals_rct_vs_nonrct_overcs.txt',sep=" ")


#Zcurve non-rcts
set.seed(2) #zcurve (Bartlos et al.) has a stochastic component and does not always converge, hence set.seed(2)
delta_zcurve_nonrct  <- zcurve_with_boot(MM_other,c=c,boots=100) 
delta_zcurve_nonrct

#Zcurve rcts
set.seed(2)
delta_zcurve_rct  <- zcurve_with_boot(MM,c=c,boots=100) 
delta_zcurve_rct

write_mm_table_paper_exact(
  out_RCT_tscores, out_other_tscores, out_RCT_articles, out_other_articles,
  floor(J_RCT_tscores), floor(J_other_tscores), floor(J_RCT_articles), floor(J_other_articles),
  eps_RCT_tscores, eps_other_tscores, eps_RCT_articles, eps_other_articles,
  t_RCT = MM$t,
  t_NE  = MM_other$t,
  tex_file = file.path(root, "output", "tables", "Table_MM_Paper.tex"),
  cv = cv,
  negate_delta = FALSE,
  # Optional: if you want to fill zcurve rows, pass these 
  zcurve_delta_RCT = delta_zcurve_rct[1],  zcurve_se_RCT = delta_zcurve_rct[2],
  zcurve_delta_NE  = delta_zcurve_nonrct[1], zcurve_se_NE = delta_zcurve_nonrct[2]
)


# full results
#table_zcurve<-cbind(delta_zcurve_rct,delta_zcurve_nonrct )
#write.csv(table_zcurve, file = 'results_zcurve.csv')

# The median number of observations:
#library(dplyr)
#median_val <- MM_data %>%
#  mutate(obs_num = as.numeric(obs)) %>%     # convert to numeric
#  distinct(title, .keep_all = TRUE) %>%     # keep one row per title
#  summarise(med = median(obs_num, na.rm = TRUE)) %>%
#  pull(med)
#median_val


#display main table results
#table_econ <- cbind(out_RCT_tscores,out_other_tscores,out_RCT_articles,out_other_articles)

#Save tables
#setwd(paste0(root,"/output/tables"))
#write.csv(table_econ, file = 'results_MM.csv')

#display results nicely
#ses<- sqrt(c(out_RCT_tscores$varest_delta,out_other_tscores$varest_delta,out_RCT_articles$varest_delta,out_other_articles$varest_delta ))
#deltas <- c(out_RCT_tscores$deltahat,out_other_tscores$deltahat,out_RCT_articles$deltahat,out_other_articles$deltahat )
#citop <- deltas+1.96*ses
#cibot <- deltas-1.96*ses
#results <- -cbind(deltas, ses,cibot,citop)
#round(results,3)