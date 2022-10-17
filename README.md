# Estimating the extrinsic incubation period of malaria using a mechanistic model of sporogony

## Authors: Isaac J Stopard, Thomas S Churcher & Ben Lambert

This GitHub repository provides all code necessary to run the multiscale model of sporogony (m:sos:) applied in "Estimating the extrinsic incubation period of malaria using a mechanistic model of sporogony". A mathematica file (*mSOS_model_mathematica.nb*) that steps through the derivations used in the model is also included. We ran the model using R version 3.6.3 and Stan version 2.21.0.

### Version 1 Release

#### Last updated: 14/09/2020

#### Stan models ("model" folder)

​	:one: *hierarchical_mSOS_all_temperature_model.stan* - implements the model fit to all the data simultaneously. Note to reduce computational time the likelihood for each unique combination of datapoint values is only calculated once.

​	:two: *mSOS_single_temperature_model_O_spz.stan* - implements the model fit to the data at 27°C, which includes the time from inoculation to oocyst and oocyst to sporozoite.

​	:three: *mSOS_single_temperature_model_spz.stan* - implements the model fit the data at the remaining temperatures (the time from inoculation to sporozoite only).

#### Running

​	:one: *run_mSOS.R* - provides R code to run the Stan models

#### Helper functions ("utils" folder)

​	:one: *data_wrangling_functions* - R helper functions to wrangle the data into the correct format for the Stan models.

​	:two: *plotting_functions.R* - R helper functions to visualise the data and model outputs.

​	:three: *plot_mSOS.R* - R script to visualise the data and model outputs.

#### R project

​	:one: *mSOS.Rproj* - R project with all necessary files and folders to run the model using original data. The "results" and "figures" folders are empty but included in the repository as locations to save the results.

### Notes

:warning: Please note that the Stan outputs are not saved on the Github repository, so run_mSOS.R will need to be run before attempting to visualise the results.

:warning: This code is released with no support and we cannot endorse any outputs from the model other than those we generate.

:warning: This model is under development, so code may change without notice.

:warning: No liability is accepted by the authors.

