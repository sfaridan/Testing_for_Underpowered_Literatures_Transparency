# Testing for Underpowered Literatures
# Stefan Faridani
# February 26, 2024
# See README for replication instructions
# This script is for transparency purposes only. 
# This code will likely change prior to publication of the article.
# Run this one!

rm(list=ls())

root <- "C:/Users/stefa/OneDrive/Documents/R/Underpowered_Literatures_Transparency"
setwd(root)

#load in functions
source("code/functions/Functions_TFUL.R")

#Generate tables and figures
source("code/scripts/Application_MM.R")
source("code/scripts/Application_ML.R")
source("code/scripts/Application_ML_Alternative.R")
source("code/scripts/compute_pvalues")