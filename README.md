---
output:
  html_document: default
  pdf_document: default
---
# Estimation and Evaluation for Generalized Quasi-Linear Mixed Model (GQLMM)

This repository contains all of the code needed to conduct the simulation studies and to reproduce the data analysis in the manuscript entitled "Quasi-Linear Predictor Function for a Novel Mixed-Effects Models".

## "code" Folder

This folder contains all of the R functions and scripts needed to implement the proposed estimation algorithm and the conditional AIC to generate the tables and figures included in the manuscript. 
The scripts require that the following packages are installed in R: *MASS*, *lme4*, *Hmisc*, and *latex2exp*. 
Once these packages are installed (the scripts will load them as needed), and once the working directory has been set to the location of the main repository on the local machine, the code should run without error. 

+ **functions.R**: Contains all functions needed to estimate the regression parameter and to calculate the conditional AIC. 
This script is sourced by all other simulation and data analaysis scripts, and so should not need to be run separately.

+ **sample-simulation-binary.R, sample-simulation-pois.R and sample-simulation-binary-misspe.R**:  Conduct simulations for the settings of the true parameter values, the cluster size and the sample sizes for each cluster.
This script generates the estimates of the regression parameters and the conditional AIC in a "result/XXX" subfolder where the output folder is automatically created.
The result files are needed to produce the figures and tables in the manuscript.
**sample-simulation-binary.R** implements the simulation for the binary response.
**sample-simulation-pois.R** implements the simulation for the Poisson response.
**sample-simulation-pois.R** conducts the simulation for the binary response reproducing the applied model misspecifies the true model, that is, omits the random-slope.

+ **figure-table-generation.R**: Imports simulation results from the "result/XXX" subfolder and produces the figures and tables in the manuscript and supporting information.
All figures will be saved as EPS format in the subfolder and all tables will be printed to the R console.

+ **data-application.R**: Conducts real data analysis.
This script generates the estimates of the regression parameters and the conditional AIC in a "result/XXX" subfolder where the output folder is automatically created.
This script also implements the backward stepwise variable selection based on the conditional AIC and the bootstrap procedure to obtain the bootstrap confidence intervals of the regression parameters and the log-odds-ratios.

## "data" Folder

This folder contains the real data imported in the script **data-application.R**.
**(All of the contents of this folder were deleted due to the copyrights.)**


## "result" Folder

This folder contains the result imported in the scripts **figure-table-generation.R**.
