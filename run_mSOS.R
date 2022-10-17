library(rstan); library(tidyverse); 

rm(list = ls())

########################################################
################ data wrangling ########################
########################################################

source("utils/data_wrangling_functions.R")

malaria_data <- read_malaria_data("data/mSOS_all_parasite_data.csv") # malaria parasite data

survival_data <- read_survival_data("data/mSOS_all_survival_data.csv") # mosquito (A. stephensi) survival data

################ all temperature model data indexing ################
placeholder <- all_temperature_data_indexing(malaria_data, survival_data)

unique_temp <- placeholder$unique_temp
unique_temp_scaled <- placeholder$unique_temp_scaled
unique_experiments <- placeholder$unique_experiments
survival_data <- placeholder$survival_data

######## oocyst data ########
oocyst_data <- subset(placeholder$malaria_data, lifestage == "cyst")

oocyst_prevalence_data <- subset(oocyst_data, ref == "Da") # Burkina Faso data - only using prevalence
oocyst_intensity_data <- subset(oocyst_data, ref != "Da")
oocyst_totals <- generate_prevalence(oocyst_prevalence_data, 1)

######## oocyst intensity data ######## 
oocyst_intensity_indices <- oocyst_intensity_indexing(oocyst_intensity_data) # getting the unique possible combinations of oocyst load, day of dissection and temperature

################ sporozoite data ################
# sporozoite prevalence data
sporozoite_data <- subset(placeholder$malaria_data, lifestage == "gland-spz")
sporozoite_totals <- generate_prevalence(sporozoite_data, 1)

################ survival data index ##################
# getting the unique possible combinations of day recorded, temperature, initial infection statuse (infected or uninfected) and status of recording (0 - censored, 1 - dead)

survival_indices <- survival_data_indexing(placeholder$survival_data)

# posterior predictive distribution times
PPD_times <- seq(0,30,0.1)
length_ppd_times <- length(PPD_times)

###############################################################
################ all temperature model ########################
###############################################################

data_in <- list(length_oocysts = as.integer(nrow(oocyst_totals)),
                oocyst_total_sampled = as.integer(oocyst_totals$sample), 
                oocyst_total_positive = as.integer(oocyst_totals$positive), 
                oocyst_index_temp = as.integer(oocyst_totals$index_temp),
                oocyst_exp = as.integer(oocyst_totals$experiment),
                oocyst_time = as.double(oocyst_totals$day_post_inf),
                
                length_sporozoites = as.integer(nrow(sporozoite_totals)), 
                sporozoite_total_sampled = as.integer(sporozoite_totals$sample),
                sporozoite_total_positive = as.integer(sporozoite_totals$positive), 
                sporozoite_index_temp = as.integer(sporozoite_totals$index_temp),
                sporozoite_exp = as.integer(sporozoite_totals$experiment),
                sporozoite_time = as.double(sporozoite_totals$day_post_inf), 
                
                length_unique_temp = as.integer(length(unique_temp_scaled)), 
                unique_temp = unique_temp_scaled,
                length_unique_exp = as.integer(nrow(unique_experiments)),
                
                length_ppd_times = as.integer(length_ppd_times), 
                PPD_times = PPD_times,
                
                length_survival_data = as.integer(nrow(survival_data)),
                survival_index = survival_indices$survival_index, 
                length_unique_survival_indices = as.integer(nrow(survival_indices$unique_survival)),
                
                unique_survival_times = as.double(survival_indices$unique_survival[,"day"]), 
                unique_survival_temp_index = as.integer(survival_indices$unique_survival[,"index_temp"]),
                unique_survival_exp = as.integer(survival_indices$unique_survival[,"experiment"]),
                unique_survival_infection_status = as.integer(survival_indices$unique_survival[,"initial_infection"]),
                unique_survival_censored = as.integer(survival_indices$unique_survival[,"survival_status"]),
                
                length_oocyst_intensity = as.integer(nrow(oocyst_intensity_data)),
                oocyst_intensity_index = as.integer(oocyst_intensity_indices$oocyst_intensity_index),
                length_unique_oocyst_intensity_indices = as.integer(nrow(oocyst_intensity_indices$unique_oocyst_intensity)),
                unique_oocyst_intensity_time = as.double(oocyst_intensity_indices$unique_oocyst_intensity[,"day_post_inf"]),
                unique_oocyst_intensity = as.double(oocyst_intensity_indices$unique_oocyst_intensity[,"number_parasites"]),
                unique_oocyst_intensity_temp_index = as.integer(oocyst_intensity_indices$unique_oocyst_intensity[,"index_temp"]),
                unique_oocyst_intensity_exp = as.integer(oocyst_intensity_indices$unique_oocyst_intensity[,"experiment"]),
                intercept = 0
)

# setting stan to run in parallel and only recompile if something has changed
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

iterations <- 4500
warmup <- 1500
chains <- 4

# # MCMC intialisation values
finit <- function(){
  list(shape_oocyst = 10, m_rate_oocyst = 0.5, c_rate_oocyst = 4, # oocyst development parameters
       shape_sporozoite = 50, rate_sporozoite = 3, # sporozoite development parameters
       mu_NB = 10, k_NB = 3, # parasite load parameters
       m_delta = 0, c_delta = 1.7, # probability of viable infection parameters
       a = 0.05, b = 0.05, beta_inf = 0.1, beta_temp = 0, #intercept = -1, # survival parameters
       error_delta = rep(0, nrow(placeholder$unique_experiments)), error_survival = rep(0, nrow(placeholder$unique_experiments)),
       sigma_error_delta = 0.005, sigma_error_survival = 0.005)
}

fit_new <- stan(file = "model/hierarchical_mSOS_all_temperature_model.stan", data = data_in, iter=iterations, chains = chains, seed=12345,
                warmup = warmup, control=list(adapt_delta=0.9, max_treedepth = 12.5, stepsize = 0.05), init = finit)

fit_wide <- stan(file = "model/hierarchical_mSOS_all_temperature_model_wider_prior.stan", data = data_in, iter=iterations, chains = chains, seed=12345,
                warmup = warmup, control=list(adapt_delta=0.9, max_treedepth = 12.5, stepsize = 0.05), init = finit)

saveRDS(fit_new, file = "results/hierarchical_mSOS_all_temperature_model_fit")

saveRDS(fit_wide, file = "results/hierarchical_mSOS_all_temperature_model_fit_wide")

###################################################################
################ single temperature models ########################
###################################################################

# MCMC intialisation values
finit_single_s <- function(){ 
  list(shape_sporozoite = 12, rate_sporozoite = 3,
       mu_NB = 10, k_NB = 3,
       delta = 0.7,
       a = 0.05, b = 0.05,
       beta_inf = 0.1) # intercept = -1
}

finit_single_27 <- function(){
  list(shape_oocyst = 10, rate_oocyst = 2,
       shape_sporozoite = 12, rate_sporozoite = 3,
       mu_NB = 10, k_NB = 3,
       delta = 0.7,
       a = 0.05, b = 0.05, beta_inf = 0.1, beta_temp = 0) # intercept = -1
}


for(i in 1:length(unique_temp)){
  data_in_single <- generate_single_temp_data_in(unique_temp[i], oocyst_totals, sporozoite_totals, 
                                                 oocyst_intensity_data, placeholder$survival_data, 
                                                 PPD_times, length_ppd_times)
  
  if(unique_temp[i] == 33){
    next # not enough time points for 33 degrees to be run as a single temperature model
  } else if(unique_temp[i] == 27){
    
    fit_27 <- stan(file = "model/mSOS_single_temperature_model_O_spz.stan", data = data_in_single, iter=iterations, chains = chains, seed=12345,
                   warmup = warmup, control=list(adapt_delta=0.9, max_treedepth = 12.5, stepsize = 0.05), init = finit_single_27)
    
    saveRDS(fit_27, file = paste0("results/mSOS_single_temperature_model_fit_new_temp_",unique_temp[i]))
    
  } else{
    
    fit_27 <- stan(file = "model/mSOS_single_temperature_model_spz.stan", data = data_in_single, iter=iterations, chains = chains, seed=12345,
                   warmup = warmup, control=list(adapt_delta=0.9, max_treedepth = 12.5, stepsize = 0.05), init = finit_single_s)
    
    saveRDS(fit_27, file = paste0("results/mSOS_single_temperature_model_fit_new_temp_",unique_temp[i]))
  }
  
}
