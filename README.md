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
    install.packages(c("doParallel", "foreach", "KFAS", "truncdist", "doSNOW", 
                        "tscount", "VGAM", "tidyverse", "dlm", "mc2d",
                        "bayesplot", "TruncatedNormal", "mvnfast", "magrittr",
                        "lubridate", "spatstat", "wesanderson", "ddst"))

## Data

In this article, we use the warpDLM methodology to analyze counts of overdose calls due to heroin and other drugs in the the city of Cincinnati.  This is derived from the full set of incident reports to the Cincinnati Fire Department, publicly available at [this link](https://data.cincinnati-oh.gov/Safety/Cincinnati-Fire-Incidents-CAD-including-EMS-ALS-BL/vnsz-a3wp).  In addition to the data set itself, data dictionaries and descriptions of coding are also provided.  This information was used to determine which calls corresponded to overdoses.

The raw data used for this project was downloaded in February 2021 and can be found [here](Data/Cincinnati_Fire_Incidents.zip).

## Application

There are two main scripts for the application. The [first file](Code/application_data_cleaning.R) is a script which inputs the raw data file and outputs the formatted count time series of drug overdoses that we are analyzing. That formatted data can be found in the Data folder as well.  The [second file](Code/application_analysis.R) has the code to run the offline Gibbs sampler and online particle filter, as well as produce the figures from the paper. The intermediate model outputs are stored [here](Outputs/ModelResults/application). Also note that this second file requires the [helper functions](Code/helper_functions.R) script, a collection of functions used in different parts of the analysis.

There are four associated figures
* [Plot of overdose series with smoothing medians overlaid](Outputs/Figures/ODCountsTS.png)
* [Plot of particle filter effective sample size across time](Outputs/Figures/PF_ESS.png)
* [Plot of particle filter time per iteration](Outputs/Figures/PF_TPL.png)
* [Plot of particle filter forecast calibration](Outputs/Figures/calibration.png)

## Simulations

The simulation forecasts are generated in four different scripts, designed to be run in a multi-core environment as they are quite computationally intensive.
* [DGLM forecasts for INGARCH data](Code/INGARCH_forecasts_dglm.R)
* [DGLM forecasts for ZIP data](Code/ZIP_forecasts_dglm.R)
* [warpDLM forecasts for INGARCH data](Code/INGARCH_forecasts_warpDLM.R)
* [warpDLM forecasts for ZIP data](Code/ZIP_forecasts_warpDLM.R)

The forecasts are stored [here](Outputs/ModelResults/simulations), and these files are used in the [simulation analysis](Code/sim_fc_analysis.R) script to compute various metrics and create figures.  In its current form, the analysis script produces figures for the INGARCH forecasts, but you can simply change one line at the beginning to use the zero-inflated Poisson forecasts instead.  

There are four associated figures
* [Box plot of rPIT calibration p-values (INGARCH)](Outputs/Figures/cal_pval_ing.png)
* [Violin plot of log score % difference over Poisson DGLM (INGARCH)](Outputs/Figures/logscore_ing.png)
* [Box plot of rPIT calibration p-values (ZIP)](Outputs/Figures/cal_pval_zip.png)
* [Violin plot of log score % difference over Poisson DGLM (ZIP)](Outputs/Figures/logscore_zip.png)
