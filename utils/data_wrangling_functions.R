# R script 
# Author: Isaac J Stopard
# Version: 0.01 
# Last updated: 12/05/2020
# Notes: functions to put data in the correct format for use in mSOS

###########################################################################
################################ FUNCTIONS ################################
###########################################################################

################ functions to read in the required data and subset the data as required ################

read_malaria_data <- function(x){
  data <- read.csv(x,header=TRUE)
  data <- subset(data, DTR == 0) # only constant temperature data in included
  data <- data[-which(data$ref == "Murdock_Thomas" & data$vector_species == "gambiae"),] # for data from the Thomas lab only A. stephensi data is included 
  return(data)
}


# selecting A. stephensi only for the Thomas lab data
read_survival_data <- function(x){
  data <- read.csv(x,header=TRUE)
  data <- subset(data, species == "stephensi") # only A. stephensi survival data is included
  return(data)
}


# function that scales the temperature, generates the temperature index and the experiment index (stratified by temperature and study)
all_temperature_data_indexing <- function(M_data, S_data){
  
  temps <- append(M_data$temp, S_data$temp)
  
  U_temp <- sort(unique(temps))
  U_temp_scaled <- sort(unique(scale(temps, center = TRUE, scale = TRUE)[,1])) # converting temperature onto the same scale as used in the other model

  for(i in 1:length(U_temp)){
    M_data[which(M_data$temp == U_temp[i]),"temp_scaled"] <- U_temp_scaled[i]
    S_data[which(S_data$temp == U_temp[i]),"temp_scaled"] <- U_temp_scaled[i]
    
    M_data[which(M_data$temp == U_temp[i]),"index_temp"] <- i
    S_data[which(S_data$temp == U_temp[i]),"index_temp"] <- i
  }
  
  U_experiments <- unique(rbind(unique(M_data[,c("index_temp", "index_ref")]), 
                                     unique(S_data[,c("index_temp", "index_ref")])))
  
  U_experiments$experiment <- seq(1, nrow(U_experiments), 1)
  
  for(i in 1:nrow(U_experiments)){
    M_data[which(M_data$index_temp == U_experiments[i,"index_temp"] &
                         M_data$index_ref == U_experiments[i,"index_ref"]), "experiment"] <- U_experiments[i,"experiment"]
    
    S_data[which(S_data$index_temp == U_experiments[i,"index_temp"] &
                          S_data$index_ref == U_experiments[i,"index_ref"]), "experiment"] <- U_experiments[i,"experiment"]
  }
  
  out <- list("unique_temp" = U_temp, "unique_temp_scaled" = U_temp_scaled, 
              "unique_experiments" = U_experiments, "malaria_data" = M_data, "survival_data" = S_data)
  
  rm(list = c("temps", "U_temp", "U_experiments", "U_temp_scaled", "M_data", "S_data"))
  
  return(out)
}

# function that generates the total number sampled and postive
# from a list of the infection statuses of individual mosquitoes
# if stratify_experiment is no this is stratified by day post infection and temperature for the posterior predictive checking and single temperature models
# if  is yes this is stratified by day post infection, temperature and experiment

generate_prevalence <- function(data, stratify_exp){
  if(stratify_exp == FALSE){
    
    totals <- unique(data[,c("day_post_inf","index_temp", "temp")])
   
    for(i in 1:nrow(totals)){
      totals[i,"sample"] <- length(which(data[,"day_post_inf"] == totals[i, "day_post_inf"]
                                       & data[,"index_temp"] == totals[i, "index_temp"]))
      
      totals[i,"positive"] <- length(which(data[,"day_post_inf"] == totals[i, "day_post_inf"] 
                                        & data[,"presence"] > 0 & data[,"index_temp"] == totals[i, "index_temp"]))
      }
    } else if(stratify_exp == TRUE){
      
      totals <- unique(data[,c("day_post_inf","index_temp", "temp", "experiment")])
      
      for(i in 1:nrow(totals)){
        totals[i,"sample"] <- length(which(data[,"day_post_inf"] == totals[i, "day_post_inf"]
                                           & data[,"index_temp"] == totals[i, "index_temp"] &
                                             data[,"experiment"] == totals[i, "experiment"]))
        
        totals[i,"positive"] <- length(which(data[,"day_post_inf"] == totals[i, "day_post_inf"] 
                                             & data[,"presence"] > 0 & data[,"index_temp"] == totals[i, "index_temp"] &
                                               data[,"experiment"] == totals[i, "experiment"]))
      }
    }
  
  rm(i)
  
  totals <- mutate(totals, prevalence = positive / sample) # prevalence
  totals <- mutate(totals, lower = prevalence - (1.96 * sqrt(prevalence * (1 - prevalence) / sample))) # 5% CI
  totals <- mutate(totals, upper = prevalence + (1.96 * sqrt(prevalence * (1 - prevalence) / sample))) # 95% CI
  totals[which(totals$lower < 0), "lower"] <- 0 # preventing the lower confidence interval being below 0
  
  return(totals)
}

### generate oocyst intensity data quantiles
generate_oocyst_intensity_summary <- function(data){
  out <- unique(data[,c("day_post_inf", "index_temp", "temp")])
  out$mean <- rep(NA, nrow(out))
  out$median <- rep(NA, nrow(out))
  out$lower <- rep(NA, nrow(out))
  out$upper <- rep(NA, nrow(out))

  for(i in 1:nrow(out)){
    placeholder <- subset(data, day_post_inf == out[i,"day_post_inf"] & temp == out[i, "temp"])
    out[i, "mean"] <- mean(placeholder[,"number_parasites"])
    out[i, "median"] <- quantile(placeholder[,"number_parasites"], 0.5)
    out[i, "lower"] <- quantile(placeholder[,"number_parasites"], 0.025)
    out[i, "upper"] <- quantile(placeholder[,"number_parasites"], 0.975)
  }
  return(out)
}


oocyst_intensity_indexing <- function(data){
  
  U_data <- unique(data[,c("day_post_inf", "number_parasites", "index_temp", "experiment")])
  
  index <- rep(NA, nrow(data))
  
  for(i in 1:nrow(U_data)){
    index[which(data$day_post_inf == U_data[i, "day_post_inf"] &
                data$number_parasites == U_data[i, "number_parasites"] &
                data$index_temp == U_data[i, "index_temp"] &
                data$experiment == U_data[i, "experiment"])] <- i
  }
  
  out <- list("unique_oocyst_intensity" = U_data, "oocyst_intensity_index" = index)
  rm(list = c("U_data", "index"))
  return(out)
}


survival_data_indexing <- function(data){
  U_data <- unique(data[c("day", "index_temp", "initial_infection", "survival_status", "experiment")])
  U_data$index <- seq(1,nrow(U_data),1)
  
  index <- rep(NA, nrow(data))
  
  for(i in 1:nrow(U_data)){
    index[which(data$day == U_data[i,"day"] & data$index_temp == U_data[i,"index_temp"] &
                data$initial_infection == U_data[i,"initial_infection"] & 
                data$survival_status == U_data[i,"survival_status"] & 
                data$experiment == U_data[i,"experiment"])] <- U_data[i,"index"]
  }
  
  out = list("unique_survival" = U_data, "survival_index" = index)
  rm(list = c("U_data", "index"))
  return(out)

}

# function that generates the data input into Stan for the single temperature data
generate_single_temp_data_in <- function(temp_, O_totals, S_totals, O_I_data, surv_data, PPD_times, length_ppd_times){
  
  # subsetting the data so only the relevant temperature is included
  surv_data <- subset(surv_data, temp == temp_)
  S_totals <- subset(S_totals, temp == temp_)
  
  # formatting the data
  U_surv <- unique(surv_data[c("day", "initial_infection", "survival_status")])
  U_surv$index <- seq(1,nrow(U_surv),1)
  surv_index <- rep(NaN, nrow(surv_data))
  
  for(i in 1:nrow(U_surv)){
    surv_index[which(surv_data$day == U_surv[i,"day"] &
                     surv_data$initial_infection == U_surv[i,"initial_infection"] & 
                           surv_data$survival_status == U_surv[i,"survival_status"])] <- U_surv[i,"index"]
  }
  
  if(temp_ == 27){
    O_totals <- subset(O_totals, temp == temp_)
    O_I_data <- subset(O_I_data, temp == temp_)
    
    # getting the unique possible combinations of oocyst load, day of dissection and temperature
    U_O_I_data <- unique(O_I_data[,c("day_post_inf", "number_parasites")])
    
    O_I_index <- rep(NaN, nrow(O_I_data))
    
    for(i in 1:nrow(U_O_I_data)){
      O_I_index[which(O_I_data$day_post_inf == U_O_I_data[i, "day_post_inf"] &
                O_I_data$number_parasites == U_O_I_data[i, "number_parasites"])] <- i
    }
    
    data_in <- list(length_oocysts = as.integer(nrow(O_totals)),
                    oocyst_total_sampled = as.integer(O_totals$sample),
                    oocyst_total_positive = as.integer(O_totals$positive),
                    oocyst_time = as.double(O_totals$day_post_inf),
                    length_oocyst_intensity = as.integer(nrow(O_I_data)),
                    oocyst_intensity_index = as.integer(O_I_index),
                    length_unique_oocyst_intensity_indices = as.integer(nrow(U_O_I_data)),
                    unique_oocyst_intensity_time = as.double(U_O_I_data$day_post_inf),
                    unique_oocyst_intensity = as.double(U_O_I_data$number_parasites),
                    length_sporozoites = as.integer(nrow(S_totals)),
                    sporozoite_total_sampled = as.integer(S_totals$sample),
                    sporozoite_total_positive = as.integer(S_totals$positive),
                    sporozoite_time = as.double(S_totals$day_post_inf),
                    length_ppd_times = as.integer(length_ppd_times),
                    PPD_times = as.double(PPD_times),
                    length_survival_data = as.integer(length(surv_index)),
                    survival_index = as.integer(surv_index),
                    length_unique_survival_indices = as.integer(nrow(U_surv)),
                    unique_survival_times = as.double(U_surv$day),
                    unique_survival_infection_status = as.integer(U_surv$initial_infection),
                    unique_survival_censored = as.integer(U_surv$survival_status),
                    intercept = 0)
  } else{
    data_in <- list(length_sporozoites = as.integer(nrow(S_totals)),
                    sporozoite_total_sampled = as.integer(S_totals$sample),
                    sporozoite_total_positive = as.integer(S_totals$positive),
                    sporozoite_time = as.double(S_totals$day_post_inf),
                    length_ppd_times = as.integer(length_ppd_times), PPD_times = PPD_times, 
                    length_survival_data = as.integer(nrow(surv_data)), 
                    survival_index = as.integer(surv_index), 
                    length_unique_survival_indices = as.integer(nrow(U_surv)),
                    unique_survival_times = as.double(U_surv$day),
                    unique_survival_infection_status = as.integer(U_surv$initial_infection),
                    unique_survival_censored = as.integer(U_surv$survival_status),
                    intercept = 0)
  }
  return(data_in)
}














