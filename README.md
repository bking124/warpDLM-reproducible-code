# WarpDLM Code

This folder contains the code necessary to reproduce the results in the paper "Warped Dynamic Linear Models for Time Series of Counts".

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
