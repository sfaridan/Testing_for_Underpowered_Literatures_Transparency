# Testing for Underpowered Literatures
# Sept 10, 2025
# Stefan Faridani
# sfaridani6@gatech.edu
# See README for replication instructions
# This script is for transparency purposes only. 
# This code will likely change prior to publication of the article.
# Run this one!

rm(list=ls())

data <- "C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/data" #change to data location
root <- "C:/Users/stefa/Onedrive/Documents/GitHub/Testing_for_Underpowered_Literatures_Transparency" #"C:/Users/sfaridani6/Documents/GitHub/Testing_for_Underpowered_Literatures_Transparency"
setwd(root)

#load in functions
source(paste0(root,"/code/functions/Functions_TFUL.R"))

#Generate tables and figures
source(paste0(root,"/code/scripts/Simulations.R"))
source(paste0(root,"/code/scripts/Application_MM.R"))
source(paste0(root,"/code/scripts/Application_ML.R"))
source(paste0(root,"/code/scripts/Application_ML_Alternative.R"))
source(paste0(root,"/code/scripts/compute_pvalues.R"))