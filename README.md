# WarpDLM Reproducibility Guide

This repository contains the code necessary to reproduce the results in the paper ["Warped Dynamic Linear Models for Time Series of Counts"](https://arxiv.org/abs/2110.14790).  There are three main folders:

* Data: Contains raw and cleaned data 
* Code: Contains all R scripts needed to run analysis
* Outputs: Contains any model outputs as well as paper figures

## Dependencies and Environments
All analysis was run in R, version 4.0.0 or higher. The bulk of the analysis was performed on a Windows laptop, with the exception of the simulations, which were run on an Rstudio Server instance hosted on a multi-core Linux server.

There are a variety of required packages, which can be installed using the following lines of code.

    #Install packages from Github
    install.packages("remotes")
    remotes::install_github("drkowal/rSTAR")
    
    #Install other 
    install.packages(c("doParallel", "foreach", "KFAS", "coda", "truncdist", 
                        "doSNOW", "tscount", "VGAM", "tidyverse", "dlm", 
                        "mc2d", "bayesplot", "TruncatedNormal", "mvnfast",
                        "magrittr", "lubridate"))

    

## Data
In this article, we use the warpDLM methodology to analyze counts of overdose calls due to heroin and other drugs in the the city of Cincinnati.  This is derived from the full set of incident reports to the Cincinnati Fire Department, publicly available at [this link](https://data.cincinnati-oh.gov/Safety/Cincinnati-Fire-Incidents-CAD-including-EMS-ALS-BL/vnsz-a3wp).  In addition to the data set itself, data dictionaries and descriptions of coding are also provided.  This information was used to determine which calls corresponded to overdoses.

The Data folder in this repository contains [a single CSV](Data/Cincinnati_Fire_Incidents.zip), which is the downloaded data set as of February 2020, used for all analysis.

## Application Results



## Simulation Results




The root folder contains the following files:
- helper_functions.R: a collection of useful functions used often in warpDLM analyses
- INGARCH_forecasts_dglm.R: the code to run the INGARCH simulation for the DGLM methods; meant to be run in a multi-core environment
- INGARCH_forecasts_warpDLM.R: the code to run the INGARCH simulation for the warpDLM methods; meant to be run in a multi-core environment
- ZIP_forecasts_dglm.R: the code to run the zero-inflated Poisson simulation for the DGLM methods; meant to be run in a multi-core environment
- ZIP_forecasts_warpDLM.R: the code to run the zero-inflated Poisson simulation for the warpDLM methods; meant to be run in a multi-core environment
- sim_fc_analysis.R: code that uses the simulation data to produce charts shown in article; default produces INGARCH sim plots, uncomment two lines of code at beginning to get ZIP plots
- application_data_cleaning.R: Takes the raw dataset and formats it into the desired time series
- Cincinnati_Fire_Incidents.csv: the raw data from the Cincinnati open data portal (downloaded Feb 2021)
- SUTSE_functions.R: contains functions to run Gibbs sampling and particle filtering for multivariate model
- application_analysis.R: reproduces the analysis for the application

There is also a subfolder ./data, which includes folders for the simulation and the application, included so that charts and analysis can be reproduced without needing to rerun lengthy simulations or redo data cleaning.
The ./data/simulation contains the output of the four simulation forecast files, and are used as input for sim_fc_analysis.R.
The ./data/application contains the time series data formatted from the original dataset, which is used for the offline and offline analyses in application_analysis.R.
