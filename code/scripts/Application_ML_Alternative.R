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
source(paste0(root,"/code/functions/Functions_TFUL_new_inference.R"))
library(ggplot2)

setwd(paste0(root,"/data"))
filename <- "ML-_Summary_Statistics.xlsx"

cs <- sqrt(1+(0:10)/10)
deltas_alt <- 0*cs
sd_deltas_alt <- 0*cs
sd_deltas_alt_uncluster <- 0*cs
for (cc in 1:length(cs))
{
    print(paste0("It ", cc, " of ", length(cs)))

    data_in = readxl::read_excel(filename,sheet="Anchoring1")
    cpower_Anchoring1<- compute_cpower(data_in$Site, data_in$`Mean (High)`, data_in$`Mean (Low)`, data_in$`SD (High)`, data_in$`SD (Low)`, data_in$`N (High)`, data_in$`N (Low)`,c=cs[cc])
    powers <- cpower_Anchoring1
  
    data_in = readxl::read_excel(filename,sheet="Anchoring2")
    cpower_Anchoring2 <- compute_cpower(data_in$Site, data_in$`Mean (High)`, data_in$`Mean (Low)`, data_in$`SD (High)`, data_in$`SD (Low)`, data_in$`N (High)`, data_in$`N (Low)`,c=cs[cc])
    powers <- rbind(powers,cpower_Anchoring2)

    data_in = readxl::read_excel(filename,sheet="Anchoring3")
    cpower_Anchoring3 <- compute_cpower(data_in$Site, data_in$`Mean (High)`, data_in$`Mean (Low)`, data_in$`SD (High)`, data_in$`SD (Low)`, data_in$`N (High)`, data_in$`N (Low)`,c=cs[cc])
    powers <- rbind(powers,cpower_Anchoring3)
    
    data_in = readxl::read_excel(filename,sheet="Anchoring4")
    cpower_Anchoring4 <- compute_cpower(data_in$Site, data_in$`Mean (High)`, data_in$`Mean (Low)`, data_in$`SD (High)`, data_in$`SD (Low)`, data_in$`N (High)`, data_in$`N (Low)`,c=cs[cc])
    powers <- rbind(powers,cpower_Anchoring4)
    
    data_in = readxl::read_excel(filename,sheet="Gambler's Fallacy")
    cpower_Gambler <- compute_cpower(data_in$Site, data_in$`Mean (Three6)`, data_in$`Mean (Two6)`, data_in$`SD (Three6)`, data_in$`SD (Two6)`, data_in$`N (Three6)`, data_in$`N (Two6)`,c=cs[cc])
    powers <- rbind(powers,cpower_Gambler)
    
    data_in = readxl::read_excel(filename,sheet="Money Priming")
    cpower_Money <- compute_cpower(data_in$Site, data_in$`Mean (Money Primg)`, data_in$`Mean (Control)`, data_in$`SD (Money Primg)`, data_in$`SD (Control)`, data_in$`N (Money Priming)`, data_in$`N (Control)`,c=cs[cc])
    powers <- rbind(powers,cpower_Money)
    
    data_in = readxl::read_excel(filename,sheet="Imagined Contact")
    cpower_Imagined <- compute_cpower(data_in$Site, data_in$`Mean (Contact)`, data_in$`Mean (Control)`, data_in$`SD (Contact)`, data_in$`SD (Control)`, data_in$`N (Contact)`, data_in$`N (Control)`,c=cs[cc])
    powers <- rbind(powers,cpower_Imagined)
    
    data_in = readxl::read_excel(filename,sheet="Sunk Costs")
    cpower_Sunk <- compute_cpower(data_in$Site, data_in$`Mean (Paid)`, data_in$`Mean (Free)`, data_in$`SD (Paid)`, data_in$`SD (Free)`, data_in$`N (Paid)`, data_in$`N (Free)`,c=cs[cc])
    powers <- rbind(powers,cpower_Sunk)
       
    data_in = readxl::read_excel(filename,sheet="Quote Attribution")
    cpower_Quote <- compute_cpower(data_in$Site, data_in$`Mean (Liked)`, data_in$`Mean (Disliked)`, data_in$`SD (Liked)`, data_in$`SD (Disliked)`, data_in$`N (Liked)`, data_in$`N (Disliked)`,c=cs[cc])
    powers <- rbind(powers,cpower_Quote)
    
    data_in = readxl::read_excel(filename,sheet="Flag Priming")
    cpower_Flag <- compute_cpower(data_in$Site, data_in$`Mean (Flag Prime)`, data_in$`Mean (Control)`, data_in$`SD (Flag Prime)`, data_in$`SD (Control)`, data_in$`N (Flag Prime)`, data_in$`N (Control)`,c=cs[cc])
    powers <- rbind(powers,cpower_Flag)
    
    data_in = readxl::read_excel(filename,sheet="Math_Art Gender")
    data_in <- data_in[c(1:22,24:28),]
    data_in$`t (unequal var)` <- as.numeric(data_in$`t (unequal var)`,c=cs[cc])
    cpower_Math <- compute_cpower(data_in$Site, data_in$`Mean (Female)`, data_in$`Mean (Male)`, data_in$`SD (Female)`, data_in$`SD (Male)`, data_in$`N (Female)`, data_in$`N (Male)`,c=cs[cc])
    powers <- rbind(powers,cpower_Math)
   
    deltas_alt[cc]  <- mean(powers[,3]-powers[,2])
    ses <- powers[,dim(powers)[2]]
    #variance of alternative estimator assuming experiments within sites are perfectly dependent
    var_cluster <- 0
    var_uncluster <- 0
    for (i in 1:length(ses)){
      for (j in 1:length(ses)){
        var_cluster <- var_cluster + ses[i]*ses[j]
        var_uncluster <- var_uncluster + ses[i]*ses[j]*(i==j)
      }
    }
    var_cluster <- var_cluster / length(ses)
    var_uncluster <- var_uncluster / length(ses)
    sd_deltas_alt[cc] <-sqrt(var_cluster)
    sd_deltas_alt_uncluster[cc] <-sqrt(var_uncluster)
    #print(cbind(cs[1:cc]^2, deltas_alt[1:cc], sd_deltas_alt[1:cc],deltas_nonrcts[1:cc]))
}

setwd(paste0(root))
makeciplot_double_gg(cs,deltas,deltas_alt,sd_deltas,sd_deltas_alt,"Unconditional","Conditional","Many Labs Replications")
ggsave("output/figures/Main_Alt_ML.pdf")

makeciplot_triple_gg(cs,deltas_rcts,deltas,deltas_alt,sd_deltas_rcts,sd_deltas,sd_deltas_alt,"RCTs","Uncond. ML","Cond. ML","Many Labs Replications")
ggsave("output/figures/RCts_Alt_ML.pdf")

makeciplot_double_gg(cs,deltas_rcts,deltas_alt,sd_deltas_rcts,sd_deltas_alt,"RCTs","Many Labs","")
ggsave("output/figures/RCts_Alt.pdf")

#Compute p-values of differences over cs
setwd(paste0(root,"/output/tables"))
pval_alt_rcts<- 2*(1-pnorm(abs(-deltas_rcts-deltas_alt)/sqrt(sd_deltas_rcts^2+sd_deltas_alt^2)))
pval_alt_rcts_uncluster<- 2*(1-pnorm(abs(-deltas_rcts-deltas_alt)/sqrt(sd_deltas_rcts^2+sd_deltas_alt_uncluster^2)))
pval_alt_non_rcts<- 2*(1-pnorm(abs(-deltas_nonrcts-deltas_alt)/sqrt(sd_deltas_nonrcts^2+sd_deltas_alt^2)))
pval_ML_non_rcts<- 2*(1-pnorm(abs(-deltas_nonrcts-deltas)/sqrt(sd_deltas_nonrcts^2+sd_deltas^2)))
pval_ML_non_rcts<- 2*(1-pnorm(abs(-deltas_nonrcts-deltas)/sqrt(sd_deltas_nonrcts^2+sd_deltas^2)))
c_sq <- cs^2
pvals_by_cs <- data.frame(c_sq,pval_alt_rcts,pval_alt_rcts_uncluster,pval_alt_non_rcts,pval_ML_non_rcts,pval_ML_non_rcts)
pvals_by_cs <- round(pvals_by_cs,4)
write.csv(pvals_by_cs,file='pvals_rct_vs_ML_overcs.csv')

#p-values when scaling tuning parameters by number of articles 
#test for equality
pval_rcts_ML_articles<- 2*(1-pnorm(abs(out_RCT_articles$deltahat-ML_by_sites$deltahat)/sqrt(out_RCT_articles$varest_delta+ML_by_sites$varest_delta  )))
pval_nonrcts_ML_articles <- 2*(1-pnorm(abs(out_RCT_tscores$deltahat-out_other_tscores$deltahat)/sqrt(out_RCT_tscores$varest_delta+out_other_tscores$varest_delta  )))
pval_rcts_alt_articles<- 2*(1-pnorm(abs(out_RCT_articles$deltahat-out_other_articles$deltahat)/sqrt(out_RCT_articles$varest_delta+out_other_articles$varest_delta  )))
pval_nonrcts_alt_articles <- 2*(1-pnorm(abs(out_RCT_tscores$deltahat-out_other_tscores$deltahat)/sqrt(out_RCT_tscores$varest_delta+out_other_tscores$varest_delta  )))
#c("P-value of RCTS=non-RCTs (by tscores): ",pval_rcts_nonrcts_tscores, ", P-value of RCTS=non-RCTs (by articles): ",pval_rcts_nonrcts_articles,  )

