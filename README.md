Replication package for "Testing for Underpowered Literatures"
Stefan Faridani
February 26, 2024

This is a replication package for transparency purposes only. A user-friendly R package for estimation and inference is forthcoming. For now, the only use-case of this package is to verify how the results in the paper were generated. This package will continue to change until the paper is published.

The paper is available at https://www.stefanfaridani.com/research

To replicate the results in the 2/26/2024 version of the paper, perform the following steps:

1. Download the file "ML-_Summary_Statistics.xlsx" from the open data posting of Klein et al. at: https://osf.io/dmf62 and read their license 

2. Download the Brodeur et al. (2020) replication package from their open data posting at: https://www.openicpsr.org/openicpsr/project/120246/version/V1/view and read their license.

3. Place the downloaded files "ML-_Summary_Statistics.xlsx" and "MM data.dta" into the data folder

4. Change the filepath called "root" in the scripts "code/scripts/Application_MM.R" and "code/scripts/Application_ML.R" to the appropriate path for your machine. 

5. Run the scripts "code/scripts/Application_MM.R" and "code/scripts/Application_ML.R" 
