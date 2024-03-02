# Testing for Underpowered Literatures
# Stefan Faridani
# February 26, 2024
# See README for replication instructions
# This script is for transparency purposes only. 
# This code will likely change prior to publication of the article.
# Run Main.R instead of this script.


#test for equality
pval_rcts_non_articles<- pval_comp(out_RCT_articles,out_other_articles) 
pval_rcts_non_tscores <- pval_comp(out_RCT_tscores,out_other_tscores)
pval_rcts_ML_articles<- pval_comp(out_RCT_articles,ML_by_sites) 
pval_rcts_ML_tscores <- pval_comp(out_RCT_tscores,ML_by_tscores)
pval_non_ML_articles<- pval_comp(out_other_articles,ML_by_sites) 
pval_non_ML_tscores <- pval_comp(out_other_tscores,ML_by_tscores)

alt_est_cluster = list(deltahat=-delta_hat_alt, varest_delta =var_cluster )
alt_est_uncluster = list(deltahat=-delta_hat_alt, varest_delta =var )
pval_rcts_alt_articles<- pval_comp(out_RCT_articles,alt_est_cluster) 
pval_rcts_alt_tscores <- pval_comp(out_RCT_tscores,alt_est_uncluster)
pval_non_alt_articles<- pval_comp(out_other_articles,alt_est_cluster) 
pval_non_alt_tscores <- pval_comp(out_other_tscores,alt_est_uncluster)

setwd(paste0(root,"/output/tables"))
write.table(c("RCTS=non-RCTs (by articles): ",pval_rcts_non_articles, ", RCTS=non-RCTs (by tscores): ",pval_rcts_non_tscores, ", RCTS=ML (by articles): ",pval_rcts_ML_articles, ", RCTS=ML (by tscores): ",pval_rcts_ML_tscores, ", non=ML (by articles): ",pval_rcts_non_articles, ", non=ML (by tscores): ",pval_non_ML_tscores, ", RCTS=ML Alt (by articles): ",pval_rcts_alt_articles, ", RCTS=alt (by tscores): ",pval_rcts_alt_tscores, ", non=alt (by articles): ",pval_non_alt_articles, ", non=alt (by tscores): ",pval_non_alt_tscores  )
            ,file='pvals.txt',sep=" ")