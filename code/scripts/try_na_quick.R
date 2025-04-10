#Set the seed (for consistent de-rounding)
set.seed(1) 

#Replace root with appropriate path
root <- "C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures"
setwd(root)

#load in functions
source(paste0(root,"/code/functions/Functions_TFUL_new_inference.R"))


noam_data <- haven::read_dta('C:/Users/stefa/OneDrive/Documents/R/Meta_Analysis_Edu/LAYS_final_database.dta')

noam_data$tscore <- noam_data$effectsize_raw/noam_data$standarderror
tscores_ready <- noam_data$tscore[!is.na(noam_data$tscore)]
papers <- as.factor(noam_data$paper[!is.na(noam_data$tscore)])
num_tscores <- length(tscores_ready)
num_papers <- length(unique(papers))

C <- 2
D <- 1e-4
sigma_Y <- 1
cv <-  1.96 #2.575
c <- sqrt(2) #sqrt(1.5)

#NA education, scale tuning parameters with number of articles
J     <- log(D*(num_papers)^(-1/3))/log(sigma_Y^2/(1+sigma_Y^2))
eps <- C*(num_papers)^(-1/3)
out_noam <- estimator(tscores_ready,J=J,cv=cv,c=c,sigma_Y=1,bandwidth=eps,studies = papers,include_pb=1)
