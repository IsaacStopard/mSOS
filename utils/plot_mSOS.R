library(rstan); library(tidyverse); library(cowplot); library(bayesplot); library(survival); library(zipfR);

rm(list = ls())

########################################################
################ data wrangling ########################
########################################################

source("utils/data_wrangling_functions.R")

source("utils/plotting_functions.R")

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

iterations <- 4500
warmup <- 1500
chains <- 4

##############################################################
################ posterior predictive means ##################
##############################################################

################ all temperature model ##################
fit_new <- readRDS(file = "results/hierarchical_mSOS_all_temperature_model_fit")

# data wrangling
fit_df <- as.data.frame(fit_new) # converting Stan fit to data frame

fit_wide <- readRDS(file = "results/hierarchical_mSOS_all_temperature_model_fit_wide")

fit_wide_df <- as.data.frame(fit_wide) # converting Stan fit to data frame

rm(list = c("fit_wide"))

# formatting the posterior predictive means
oocyst_prop_ppd_df <- prop_ppd_function(fit_df, unique_temp, length_ppd_times, iterations, warmup, chains, length(unique_temp), "oocyst_prevalence_ppd")
sporozoite_prop_ppd_df <- prop_ppd_function(fit_df, unique_temp, length_ppd_times, iterations, warmup, chains, length(unique_temp), "sporozoite_prevalence_ppd")
oocyst_intensity_ppd_df <- prop_ppd_function(fit_df, unique_temp, length_ppd_times, iterations, warmup, chains, length(unique_temp), "oocyst_intensity_ppd")

oocyst_prop_ppd_df_wide <- prop_ppd_function(fit_wide_df, unique_temp, length_ppd_times, iterations, warmup, chains, length(unique_temp), "oocyst_prevalence_ppd")
sporozoite_prop_ppd_df_wide <- prop_ppd_function(fit_wide_df, unique_temp, length_ppd_times, iterations, warmup, chains, length(unique_temp), "sporozoite_prevalence_ppd")
oocyst_intensity_ppd_df_wide <- prop_ppd_function(fit_wide_df, unique_temp, length_ppd_times, iterations, warmup, chains, length(unique_temp), "oocyst_intensity_ppd")


# oocyst intensity data wrangling
k <- mean(rstan::extract(fit_new, pars = "k_NB")[[1]])
oocyst_intensity_ppd_df$lower_m <- qnbinom(0.05, mu = oocyst_intensity_ppd_df$mean, size = k)
oocyst_intensity_ppd_df$upper_m <- qnbinom(0.95, mu = oocyst_intensity_ppd_df$mean, size = k)

# survival data
St_infected_prop_ppd_df <- prop_ppd_function(fit_df, unique_temp, length_ppd_times, iterations, warmup, chains, length(unique_temp), "St_infected_blood_fed_ppd")
St_uninfected_prop_ppd_df <- prop_ppd_function(fit_df, unique_temp, length_ppd_times, iterations, warmup, chains, length(unique_temp), "St_uninfected_ppd")

St_uninfected_prop_ppd_df$infection_status <- rep(0, nrow(St_uninfected_prop_ppd_df))
St_infected_prop_ppd_df$infection_status <- rep(1, nrow(St_infected_prop_ppd_df))
St_ppd_df <- rbind(St_uninfected_prop_ppd_df, St_infected_prop_ppd_df) # merging so in correct format for plotting
St_ppd_df$strata_infection <- St_ppd_df$infection_status

# theme for plotting the data
blank_theme<-theme(axis.title=element_text(colour="black", size = 20),axis.text=element_text(colour="black", size = 18),
                   axis.ticks=element_line(colour="black"),panel.background = element_blank(),
                   axis.line = element_line(colour = "black"), legend.title = element_text(size = 20),
                   legend.text=element_text(size=20), legend.key = element_rect(fill = "white"))

################ traceplot ##################
parameters <- c("shape_oocyst", "m_rate_oocyst", "c_rate_oocyst",
                "shape_sporozoite", "rate_sporozoite", "mu_NB", "k_NB", 
                "m_delta", "c_delta", "a", "b", "beta_inf", "beta_temp", # "intercept", 
                "sigma_error_delta", "sigma_error_survival")

fit_bayesplot <- as.array(fit_new, pars = parameters)
parameter_labels <- c("alpha[GO]", "m[beta]", "c[beta]", "alpha[OS]",
                      "beta[OS]", "mu", "k", "m[delta]", "c[delta]", "a[s]", "b[s]",
                      "beta[E]", "beta[C]", #"C[cox]", 
                      "sigma[delta]", "sigma[survival]")

dimnames(fit_bayesplot)[[3]] <- parameter_labels

color_scheme_set("darkgray")
png(file = "figures/traceplot.png", width = 1600, height = 800)
mcmc_trace(fit_bayesplot, parameter_labels, n_warmup = 0, 
           facet_args = list(labeller = label_parsed)) + facet_text(face = "bold", size = 20) +
  theme(axis.title=element_text(colour="black", size = 18),axis.text=element_text(colour="black", size = 18),
        axis.ticks=element_line(colour="black"),panel.background = element_blank(),
        axis.line = element_line(colour = "black"), legend.title = element_text(size = 20),
        legend.text=element_text(size=20), legend.key = element_rect(fill = "white")) + 
  guides(colour = guide_legend(override.aes = list(size = 7)))
dev.off()

# pairs plots of the gamma distribution mean and variance
shape_o <-  rstan::extract(fit_new, "shape_oocyst")[[1]]
shape_s <- rstan::extract(fit_new, "shape_sporozoite")[[1]]
rate_s <- rstan::extract(fit_new, "rate_sporozoite")[[1]]

mu_total_s <- data.frame(rstan::extract(fit_new, "mu_total_sporozoite")[[1]])
colnames(mu_total_s) <- unique_temp
mu_total_s <- mu_total_s %>% gather(key = "temp", value = "mu_total_s", 1:ncol(mu_total_s))

sigma_sq_s <- data.frame(rstan::extract(fit_new, "sigma_sq_sporozoite")[[1]])
colnames(sigma_sq_s) <- unique_temp
sigma_sq_s <- sigma_sq_s %>% gather(key = "temp", value = "sigma_sq_s", 1:ncol(sigma_sq_s))

ggplot(data = cbind(mu_total_s, sigma_sq_s), aes(x = mu_total_s, y = sigma_sq_s)) +
  geom_point(shape = 21, alpha = 0.1, fill = "skyblue", col = "skyblue1") + theme_bw() +
  facet_wrap(~ temp, scales = "free")

rate_o <- rstan::extract(fit_new, "rate_oocyst")[[1]]

gamma_o_df <- bind_rows(lapply(seq(1, ncol(rate_o)), function(i, unique_temp, shape_o, rate_o){
  return(data.frame("mean" = shape_o/rate_o[,i],
                    "var" = shape_o/rate_o[,i]^2,
                    "temp" = unique_temp[i]))
  },
  unique_temp = unique_temp,
  shape_o = shape_o,
  rate_o = rate_o))

png(file = "figures/mean_var_pairs.png", height = 1100, width = 1000)
plot_grid(
  ggplot(data = gamma_o_df %>% mutate(temp = paste0(temp,"°C")), aes(x = mean, y = var)) +
    geom_point(shape = 21, alpha = 0.17, fill = "deepskyblue", col = "deepskyblue3") +
    facet_wrap(~temp, scales = "free") + blank_theme  +
    xlab("Mean of the gamma distributed\nparasite development times from G to O") +
    ylab("Variance of the gamma distributed\nparasite development times from G to O") +
    theme(strip.text.x = element_text(size = 15)),
  
  plot_grid(
    ggplot(data = data.frame("mean_s" = shape_s/rate_s, "var_s" = shape_s/rate_s^2),
       aes(x = mean_s, y = var_s)) +
  geom_point(shape = 21, col = "deepskyblue3", alpha = 0.17, fill = "deepskyblue") + blank_theme +
  scale_x_continuous(breaks = seq(7.5, 9, 0.5), limits = c(7.5, 9)) +
  xlab("Mean of the gamma distributed\nparasite development times from O to S") +
  ylab("Variance of the gamma distributed\nparasite development times from O to S"),
  NULL, nrow = 1, rel_widths = c(0.5, 0.25)
  ),
  nrow = 2, rel_heights = c(1, 0.75),
  labels = c("A", "B")
)
dev.off()

################## all temperature model fits ###################
# oocyst plot data
oocyst_totals_all <- generate_prevalence(oocyst_prevalence_data, 0)
sporozoite_totals_all <- generate_prevalence(sporozoite_data, 0)
oocyst_intensity_summary <- generate_oocyst_intensity_summary(oocyst_intensity_data)

# oocyst plots
x_axis <- ggdraw() + draw_label("Days post blood feed", color = "black", fontface = "plain", size = 20,
                                vjust = 0.001) + theme(plot.margin = margin(0, 0, 0, 0))

p <- plot_grid(oocyst_prevalence_plot(oocyst_prop_ppd_df, oocyst_totals_all, 27, unique_temp, blank_theme) + 
            theme(axis.title.x = element_blank()), 
          oocyst_intensity_plot_PR(oocyst_intensity_ppd_df, oocyst_intensity_summary, 27, blank_theme) + 
            theme(axis.title.x = element_blank()), 
          oocyst_intensity_plot_PR(oocyst_intensity_ppd_df, oocyst_intensity_summary, 30, blank_theme) +
            theme(axis.title.x = element_blank()), 
          oocyst_intensity_plot_PR(oocyst_intensity_ppd_df, oocyst_intensity_summary, 33, blank_theme) + 
            theme(axis.title.x = element_blank()),
          nrow = 2, ncol = 2)

png(file = "figures/oocysts_all_temp_model.png", width = 800, height = 500)
plot_grid(p, x_axis, ncol = 1, rel_heights = c(1, 0.05))
dev.off()

# sporozoite plots
png(file = "figures/sporozoite_prevalence_all_temp_model.png", width = 800, height = 600)
sporozoite_prevalence_plot_all(sporozoite_prop_ppd_df, sporozoite_totals_all, blank_theme)
dev.off()


sporozoite_prevalence_plot_all <- function(sporozoite_prop_ppd_df, sporozoite_prop_ppd_df_wide, sporozoite_totals_all, blank_theme){
  sporozoite_prop_ppd_df$temp_label <- paste0(sporozoite_prop_ppd_df$temp, "°C")
  sporozoite_prop_ppd_df$prior <- rep("2.5", nrow(sporozoite_prop_ppd_df))
  sporozoite_prop_ppd_df_wide$temp_label <- paste0(sporozoite_prop_ppd_df_wide$temp, "°C")
  sporozoite_prop_ppd_df_wide$prior <- rep("4.5", nrow(sporozoite_prop_ppd_df_wide))
  sporozoite_prop_ppd_df <- rbind(sporozoite_prop_ppd_df, sporozoite_prop_ppd_df_wide)
  
  sporozoite_totals_all$temp_label <- paste0(sporozoite_totals_all$temp, "°C")
  ggplot() +
    geom_ribbon(data = sporozoite_prop_ppd_df, 
                aes(x = day_post_inf, ymin = lower, ymax = upper, fill = prior), alpha = 0.3) +
    geom_line(data = sporozoite_prop_ppd_df, 
              aes(x = day_post_inf, y = median, colour = prior), linetype = 1, size = 0.75, alpha = 0.75) +
    geom_pointrange(data = sporozoite_totals_all, 
                    aes(x = day_post_inf, y = prevalence, ymin = lower, ymax = upper), col = "black") +
    facet_wrap(vars(temp_label), scales = "free") +
    xlab("Days post blood-feed") + ylab("Sporozoite prevalence") + 
    scale_y_continuous(limits = c(0,0.95), breaks = seq(0,0.8,0.2), labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(limits = c(0,30), breaks = seq(0,30,10)) + blank_theme + 
    theme(strip.text = element_text(face="bold", size=20))
}

oocyst_prevalence_plot <- function(oocyst_prop_ppd_df, oocyst_proportion_count_data, temp_, unique_temp_original, blank_theme){
  oocyst_prop_ppd_df$temp_label <- paste0(oocyst_prop_ppd_df$temp, "°C")
  oocyst_prop_ppd_df$prior <- rep("2.5", nrow(oocyst_prop_ppd_df))
  oocyst_prop_ppd_df_wide$temp_label <- paste0(oocyst_prop_ppd_df_wide$temp, "°C")
  oocyst_prop_ppd_df_wide$prior <- rep("4.5", nrow(oocyst_prop_ppd_df_wide))
  oocyst_prop_ppd_df <- rbind(oocyst_prop_ppd_df, oocyst_prop_ppd_df_wide)
  
  oocyst_proportion_count_data$temp_label <- paste0(oocyst_proportion_count_data$temp,"°C")
  
  temp_ <- 27
  
  ggplot() +
    geom_ribbon(data = subset(oocyst_prop_ppd_df, temp == temp_ & day_post_inf <= 10), 
                aes(x = day_post_inf, ymin = lower, ymax = upper, fill = prior), alpha = 0.4) +
    geom_line(data = subset(oocyst_prop_ppd_df, temp == temp_ & day_post_inf <= 10), 
              aes(x = day_post_inf, y = median, colour = prior), linetype = 1, size = 0.75) +
    geom_pointrange(data = subset(oocyst_proportion_count_data, index_temp == which(unique_temp_original == temp_)), 
                    aes(x = day_post_inf, y = prevalence, ymin = lower, ymax = upper), colour = "black") +
    xlab("Days post blood-feed") + ylab("Oocyst prevalence") + 
    scale_y_continuous(limits = c(0,0.95), breaks = seq(0,0.8,0.2), labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(limits = c(0,10), breaks = seq(0,10,2)) + blank_theme +
    facet_grid(. ~ temp_label) +  theme(strip.text = element_text(face="bold", size=20))
}


########################################
######### logistic growth model ########
########################################
day_cut_off <- rep(NA, length(unique_temp))
for(i in 1:length(unique_temp)){
  placeholder_df <- subset(sporozoite_totals_all, temp == unique_temp[i])
  day_cut_off[i] <- placeholder_df[which(placeholder_df$prevalence == max(placeholder_df$prevalence)),"day_post_inf"]
}

l_fit_21 <- nls(presence ~ g/(1+exp(-k*(day_post_inf-t))),data=subset(sporozoite_data, temp == unique_temp[1] & day_post_inf <= day_cut_off[1]),start=list(g=0.4,k=2,t=12), algorithm = "port")
l_fit_24 <- nls(presence ~ g/(1+exp(-k*(day_post_inf-t))),data=subset(sporozoite_data, temp == unique_temp[2] & day_post_inf <= day_cut_off[2]),start=list(g=0.4,k=2,t=12), algorithm = "port")
l_fit_27 <- nls(presence ~ g/(1+exp(-k*(day_post_inf-t))),data=subset(sporozoite_data, temp == unique_temp[3] & day_post_inf <= day_cut_off[3]),start=list(g=0.4,k=2,t=12), algorithm = "port")
l_fit_30 <- nls(presence ~ g/(1+exp(-k*(day_post_inf-t))),data=subset(sporozoite_data, temp == unique_temp[4] & day_post_inf <= day_cut_off[4]),start=list(g=0.7,k=0.25,t=10), algorithm = "port")
l_fit_32 <- nls(presence ~ g/(1+exp(-k*(day_post_inf-t))),data=subset(sporozoite_data, temp == unique_temp[5] & day_post_inf <= day_cut_off[5]),start=list(g=0.51,k=0.44,t=7.8), algorithm = "port")
l_fit_34 <- nls(presence ~ g/(1+exp(-k*(day_post_inf-t))),data=subset(sporozoite_data, temp == unique_temp[7] & day_post_inf <= day_cut_off[7]), start=list(g=0.2,k=1,t=7), algorithm = "port")

logistic_fit_df <- rbind(data.frame("temp" = rep(21, length(PPD_times)), "day_post_inf" = PPD_times, 
                                    "median" = summary(l_fit_21)$coefficients["g",1] / (1 + exp(- summary(l_fit_21)$coefficients["k",1] *
                                                                                                    (PPD_times - summary(l_fit_21)$coefficients["t",1])))),
                         data.frame("temp" = rep(24, length(PPD_times)), "day_post_inf" = PPD_times, 
                                    "median" = summary(l_fit_24)$coefficients["g",1] / (1 + exp(- summary(l_fit_24)$coefficients["k",1] *
                                                                                                  (PPD_times - summary(l_fit_24)$coefficients["t",1])))),
                         data.frame("temp" = rep(27, length(PPD_times)), "day_post_inf" = PPD_times, 
                                    "median" = summary(l_fit_27)$coefficients["g",1] / (1 + exp(- summary(l_fit_27)$coefficients["k",1] *
                                                                                                  (PPD_times - summary(l_fit_27)$coefficients["t",1])))),
                         data.frame("temp" = rep(30, length(PPD_times)), "day_post_inf" = PPD_times, 
                                    "median" = summary(l_fit_30)$coefficients["g",1] / (1 + exp(- summary(l_fit_30)$coefficients["k",1] *
                                                                                                  (PPD_times - summary(l_fit_30)$coefficients["t",1])))),
                         data.frame("temp" = rep(32, length(PPD_times)), "day_post_inf" = PPD_times, 
                                    "median" = summary(l_fit_32)$coefficients["g",1] / (1 + exp(- summary(l_fit_32)$coefficients["k",1] *
                                                                                                  (PPD_times - summary(l_fit_32)$coefficients["t",1])))),
                         data.frame("temp" = rep(34, length(PPD_times)), "day_post_inf" = PPD_times, 
                                    "median" = summary(l_fit_34)$coefficients["g",1] / (1 + exp(- summary(l_fit_34)$coefficients["k",1] *
                                                                                                  (PPD_times - summary(l_fit_34)$coefficients["t",1])))))
png(file = "figures/logistic_all_temp_model.png", width = 1400, height = 600)
logistic_fit_plot_all(sporozoite_totals_all, logistic_fit_df, unique_temp,
                      day_cut_off, blank_theme)
dev.off()

# survival data
# Kaplan 
surv_frame_21 <- createSurvivalFrame(survfit(Surv(day, survival_status) ~ initial_infection, data = subset(survival_data, temp == 21)))
surv_frame_21$temp <- rep(21, nrow(surv_frame_21))
surv_frame_24 <- createSurvivalFrame(survfit(Surv(day, survival_status) ~ initial_infection, data = subset(survival_data, temp == 24)))
surv_frame_24$temp <- rep(24, nrow(surv_frame_24))
surv_frame_27 <- createSurvivalFrame(survfit(Surv(day, survival_status) ~ initial_infection, data = subset(survival_data, temp == 27)))
surv_frame_27$temp <- rep(27, nrow(surv_frame_27))
surv_frame_30 <- createSurvivalFrame(survfit(Surv(day, survival_status) ~ initial_infection, data = subset(survival_data, temp == 30)))
surv_frame_30$temp <- rep(30, nrow(surv_frame_30))
surv_frame_32 <- createSurvivalFrame(survfit(Surv(day, survival_status) ~ initial_infection, data = subset(survival_data, temp == 32)))
surv_frame_32$temp <- rep(32, nrow(surv_frame_32))
surv_frame_33 <- createSurvivalFrame(survfit(Surv(day, survival_status) ~ initial_infection, data = subset(survival_data, temp == 33)))
surv_frame_33$temp <- rep(33, nrow(surv_frame_33))
surv_frame_34 <- createSurvivalFrame(survfit(Surv(day, survival_status) ~ initial_infection, data = subset(survival_data, temp == 34)))
surv_frame_34$temp <- rep(34, nrow(surv_frame_34))
surv_frame <- rbind(surv_frame_21, surv_frame_24, surv_frame_27, surv_frame_30, surv_frame_32, surv_frame_33, surv_frame_34)

png(file = "figures/survival_all_temp_model.png", width = 1000, height = 600)
survival_plot(surv_frame, St_ppd_df, blank_theme)
dev.off()

png(file = "figures/kaplan_meier_plot.png", width = 1000, height = 600)
kaplan_meier_plot(surv_frame, blank_theme)
dev.off()

################ single temperature models ##################
indexes <- c(1,2,3,4,5,7) # indexes to exclude 33

# reading in the RDS model fit files
for(i in 1:length(indexes)){
  index <- indexes[i]
  assign(paste0("fit_",unique_temp[index]), 
         readRDS(file = paste0("results/mSOS_single_temperature_model_fit_new_temp_",unique_temp[index])))
}


oocyst_prop_ppd_df_27 <- prop_ppd_function_single(length_ppd_times, iterations, warmup, chains, fit_27, "oocyst_prevalence_ppd")
oocyst_intensity_ppd_df_27 <- prop_ppd_function_single(length_ppd_times, iterations, warmup, chains, fit_27, "oocyst_intensity_ppd")

# calculating the PPD for intensity
k_27 <- mean(rstan::extract(fit_27, par = "k_NB")[[1]])
oocyst_intensity_ppd_df_27$lower_m <- qnbinom(0.05, mu = oocyst_intensity_ppd_df_27$mean, size = k_27)
oocyst_intensity_ppd_df_27$upper_m <- qnbinom(0.95, mu = oocyst_intensity_ppd_df_27$mean, size = k_27)

for(i in 1:length(indexes)){
  assign(paste0("sporozoite_prop_ppd_df_",unique_temp[indexes[i]]), 
         prop_ppd_function_single(length_ppd_times, iterations, warmup, chains, get(paste0("fit_", unique_temp[indexes[i]])), "sporozoite_prevalence_ppd"))
  
  assign(paste0("St_infected_prop_ppd_df_", unique_temp[indexes[i]]),
         prop_ppd_function_single(length_ppd_times, iterations, warmup, chains, get(paste0("fit_", unique_temp[indexes[i]])), "St_infected_blood_fed_ppd"))
  
  assign(paste0("St_uninfected_prop_ppd_df_", unique_temp[indexes[i]]),
         prop_ppd_function_single(length_ppd_times, iterations, warmup, chains, get(paste0("fit_", unique_temp[indexes[i]])), "St_uninfected_ppd"))
}

St_uninfected_prop_ppd_df_21$infection_status <- rep(0, nrow(St_uninfected_prop_ppd_df_21))
St_uninfected_prop_ppd_df_24$infection_status <- rep(0, nrow(St_uninfected_prop_ppd_df_24))
St_uninfected_prop_ppd_df_27$infection_status <- rep(0, nrow(St_uninfected_prop_ppd_df_27))
St_uninfected_prop_ppd_df_30$infection_status <- rep(0, nrow(St_uninfected_prop_ppd_df_30))
St_uninfected_prop_ppd_df_32$infection_status <- rep(0, nrow(St_uninfected_prop_ppd_df_32))
St_uninfected_prop_ppd_df_34$infection_status <- rep(0, nrow(St_uninfected_prop_ppd_df_34))

St_infected_prop_ppd_df_21$infection_status <- rep(1, nrow(St_infected_prop_ppd_df_21))
St_infected_prop_ppd_df_24$infection_status <- rep(1, nrow(St_infected_prop_ppd_df_24))
St_infected_prop_ppd_df_27$infection_status <- rep(1, nrow(St_infected_prop_ppd_df_27))
St_infected_prop_ppd_df_30$infection_status <- rep(1, nrow(St_infected_prop_ppd_df_30))
St_infected_prop_ppd_df_32$infection_status <- rep(1, nrow(St_infected_prop_ppd_df_32))
St_infected_prop_ppd_df_34$infection_status <- rep(1, nrow(St_infected_prop_ppd_df_34))

St_prop_ppd_df_21 <- rbind(St_uninfected_prop_ppd_df_21, St_infected_prop_ppd_df_21)
St_prop_ppd_df_21$temp <- rep(21, nrow(St_prop_ppd_df_21))
St_prop_ppd_df_24 <- rbind(St_uninfected_prop_ppd_df_24, St_infected_prop_ppd_df_24)
St_prop_ppd_df_24$temp <- rep(24, nrow(St_prop_ppd_df_24))
St_prop_ppd_df_27 <- rbind(St_uninfected_prop_ppd_df_27, St_infected_prop_ppd_df_27)
St_prop_ppd_df_27$temp <- rep(27, nrow(St_prop_ppd_df_27))
St_prop_ppd_df_30 <- rbind(St_uninfected_prop_ppd_df_30, St_infected_prop_ppd_df_30)
St_prop_ppd_df_30$temp <- rep(30, nrow(St_prop_ppd_df_30))
St_prop_ppd_df_32 <- rbind(St_uninfected_prop_ppd_df_32, St_infected_prop_ppd_df_32)
St_prop_ppd_df_32$temp <- rep(32, nrow(St_prop_ppd_df_32))
St_prop_ppd_df_34 <- rbind(St_uninfected_prop_ppd_df_34, St_infected_prop_ppd_df_34)
St_prop_ppd_df_34$temp <- rep(34, nrow(St_prop_ppd_df_34))

St_prop_ppd_df_single_temp_model <- rbind(St_prop_ppd_df_21, St_prop_ppd_df_24, St_prop_ppd_df_27,
                                          St_prop_ppd_df_30, St_prop_ppd_df_32, St_prop_ppd_df_34)
surv_frame_single_temp_model <- surv_frame[-which(surv_frame$temp == 33),]

# survival plots
png(file = "figures/survival_single_temp_model.png", width = 1100, height = 450)
survival_plot(surv_frame_single_temp_model, St_prop_ppd_df_single_temp_model, blank_theme)
dev.off()

# 27 plots
x_axis <- ggdraw() + draw_label("Days post blood feed", color = "black", fontface = "plain", size = 20,
                                vjust = 0.001) + theme(plot.margin = margin(0, 0, 0, 0))
p <- plot_grid(oocyst_prevalence_plot_single(oocyst_prop_ppd_df_27, oocyst_totals_all, 27, blank_theme, "A") + 
                 theme(axis.title.x = element_blank()),
               oocyst_intensity_plot_single_PR(oocyst_intensity_ppd_df_27, oocyst_intensity_summary, 27, blank_theme, "B") +
                 theme(axis.title.x = element_blank()),
               sporozoite_prevalence_plot_single(sporozoite_prop_ppd_df_27, sporozoite_totals_all, 27, blank_theme, "C") + 
                 theme(axis.title.x = element_blank()),
               nrow = 2, ncol = 2)

png(file = "figures/oocysts_single_temp_model_27.png", width = 800, height = 500)
plot_grid(p, x_axis, ncol = 1, rel_heights = c(1, 0.05))
dev.off()

###################################################################################
########## getting the values for the results for the 27 degrees celsius ##########
###################################################################################

# checking the values
# times are only calculated to the nearest 0.1 days
# the model outputs are rounded to the 3 decimal places and the time-point nearest the desired model output is selected
oocyst_check <- round(oocyst_prop_ppd_df_27$upper, digits = 4)
oocyst_prop_ppd_df_27[which(abs(oocyst_check - 0.1) == min(abs(oocyst_check - 0.1))),"day_post_inf"]
oocyst_prop_ppd_df_27[which(abs(oocyst_check - max(oocyst_check)) == min(abs(oocyst_check - max(oocyst_check)))),"day_post_inf"]
max(oocyst_check)

oocyst_intensity_check <- round(oocyst_intensity_ppd_df_27$upper, digits = 4)
oocyst_intensity_ppd_df_27[which(abs(oocyst_intensity_check - max(oocyst_intensity_check)) == min(abs(oocyst_intensity_check - max(oocyst_intensity_check)))),"day_post_inf"]
max(oocyst_intensity_check)

sporozoite_check <- round(sporozoite_prop_ppd_df_27$upper, digits = 4)
a <- subset(sporozoite_prop_ppd_df_27, day_post_inf < sporozoite_prop_ppd_df_27[max(which(sporozoite_check == max(sporozoite_check))),"day_post_inf"])
a_check <- round(a$median, digits = 4)
sporozoite_prop_ppd_df_27[which(abs(sporozoite_check - 0.1) == min(abs(sporozoite_check - 0.1))),"day_post_inf"]
a[which(abs(a_check - max(sporozoite_check)) == min(abs(a_check - max(sporozoite_check)))),"day_post_inf"]
max(sporozoite_check)

############################################################################################
####################### sporozoite plots - single temperature models #######################
############################################################################################

sporozoite_prop_ppd_df_21$temp <- rep(21, nrow(sporozoite_prop_ppd_df_21))
sporozoite_prop_ppd_df_24$temp <- rep(24, nrow(sporozoite_prop_ppd_df_24))
sporozoite_prop_ppd_df_27$temp <- rep(27, nrow(sporozoite_prop_ppd_df_27))
sporozoite_prop_ppd_df_30$temp <- rep(30, nrow(sporozoite_prop_ppd_df_30))
sporozoite_prop_ppd_df_32$temp <- rep(32, nrow(sporozoite_prop_ppd_df_32))
sporozoite_prop_ppd_df_34$temp <- rep(34, nrow(sporozoite_prop_ppd_df_34))

sporozoite_prop_ppd_df_single <- rbind(sporozoite_prop_ppd_df_21, sporozoite_prop_ppd_df_24,
                                       sporozoite_prop_ppd_df_27, sporozoite_prop_ppd_df_30, 
                                       sporozoite_prop_ppd_df_32, sporozoite_prop_ppd_df_34)

sporozoite_totals_single <- sporozoite_totals_all[-which(sporozoite_totals_all$temp == 33),]

png(file = "figures/sporozoite_prevalence_single_temp_model.png", width = 1100, height = 450)
sporozoite_prevalence_plot_all(sporozoite_prop_ppd_df_single, sporozoite_totals_single, blank_theme)
dev.off()

############################################
################ EIP plot ##################
############################################

# creating a data frame of the parameter values for each iteration
fit_extract <- rstan::extract(fit_new)

parameter_df <- data.frame(iteration = seq(1:(iterations - warmup) * chains),
                           shape_oocyst = fit_extract$shape_oocyst,
                           m_rate_oocyst = fit_extract$m_rate_oocyst,
                           c_rate_oocyst = fit_extract$c_rate_oocyst,
                           shape_sporozoite = fit_extract$shape_sporozoite,
                           rate_sporozoite = fit_extract$rate_sporozoite,
                           mu = fit_extract$mu_NB,
                           k = fit_extract$k_NB)
# required inputs
time_input <- seq(0.0001, 50, 0.0001)
temp_input <- round(seq(unique_temp_scaled[1], unique_temp_scaled[7], (unique_temp_scaled[7] - unique_temp_scaled[6])/2), digits = 9)
unique_temp_scaled <- round(unique_temp_scaled, digits = 9)
iteration_index <- seq(5,((iterations - warmup) * chains), 5) # thinning the iterations
percentiles <- c(0.1, 0.5, 0.9) # EIP percentile
quantiles <- c(0.05, 0.5, 0.95)

# estimating the EIP for every 5th iteration - takes a long time to run
to_go <- generate_EIP_values(parameter_df, time_input, temp_input, unique_temp_scaled, iteration_index, percentiles) # need to run this again is the RDS hasn't been saved
to_go <- readRDS(file = "results/hierarchical_mSOS_EIP_values")
to_go_percentiles_plot <- generate_EIP_quantiles(to_go, temp_input, percentiles, quantiles) 
iteration_index_single <- seq(1,((iterations - warmup)*chains), 1)
to_go_single <- generate_EIP_values_single(unique_temp, iteration_index_single, time_input, percentiles) # need to run this again is the RDS hasn't been saved
to_go_single_all <- readRDS(file = "results/mSOS_EIP_values_single")

single_all_range_data <- expand.grid(unique(to_go_single_all$percentile), unique(to_go_single_all$temps_single_all))
colnames(single_all_range_data) <- c("EIP_X", "temp")
single_all_range_data$lower = rep(NaN, nrow(single_all_range_data))
single_all_range_data$median = rep(NaN, nrow(single_all_range_data))
single_all_range_data$upper = rep(NaN, nrow(single_all_range_data))

for(i in 1:nrow(single_all_range_data)){
  a <- subset(to_go_single_all, percentile == single_all_range_data[i,"EIP_X"] 
              & temps_single_all == single_all_range_data[i,"temp"])
  single_all_range_data[i,"lower"] <- quantile(a$EIP_times_single_all, 0.05)[[1]]
  single_all_range_data[i,"median"] <- quantile(a$EIP_times_single_all, 0.5)[[1]]
  single_all_range_data[i,"upper"] <- quantile(a$EIP_times_single_all, 0.95)[[1]]
}

for(i in 1:length(unique_temp)){
  single_all_range_data[which(single_all_range_data$temp == unique_temp[i]), "temp"] <- unique_temp_scaled[i]
}

# probability of viable infection (realised transmission probability)
to_go_delta_percentiles_plot <- extract_delta(fit_extract, iteration_index_single, temp_input, quantiles)
saveRDS(to_go_delta_percentiles_plot, file = "results/hierarchical_mSOS_delta_values")
to_go_delta_percentiles_plot <- readRDS(file = "results/hierarchical_mSOS_delta_values")

delta_single_values <- extract_delta_single(indexes, iteration_index_single, unique_temp_scaled, unique_temp)

EIP_plot <- EIP_percentile_plot(to_go_percentiles_plot, single_all_range_data)

delta_plot_ <- delta_plot(to_go_delta_percentiles_plot, delta_single_values)

x_axis <- ggdraw() + draw_label("Temperature / °C", color = "black", fontface = "plain", size = 20,
                                vjust = 0.01) + theme(plot.margin = margin(0, 0, 0, 0))
p <- plot_grid(EIP_plot + theme(axis.title.x = element_blank(), legend.key.size = unit(1.25,"cm")) + 
                                  guides(shape = guide_legend(override.aes = list(size = 1))), 
               delta_plot_ + theme(axis.title.x = element_blank()))


png(file = "figures/EIP_delta_plot.png", width = 1100, height = 450)
plot_grid(p, x_axis, ncol = 1, rel_heights = c(1, 0.05))
dev.off()

# logistic fits EIP estimates
l_EIP <- single_all_range_data
l_EIP$lower <- rep(0, nrow(l_EIP))
l_EIP$median <- rep(NA, nrow(l_EIP))
l_EIP$upper <- rep(0, nrow(l_EIP))
for(i in 1:nrow(l_EIP)){
  a <- summary(get(paste0("l_fit_",unique_temp[which(l_EIP[i,"temp"] == unique_temp_scaled)])))
  
  l_EIP[i, "median"] <- -((log((a$coefficients["g",1] / (a$coefficients["g",1] * l_EIP[i, "EIP_X"])) - 1)) 
    / a$coefficients["k",1]) + a$coefficients["t",1]
}

l_delta <- delta_single_values
l_delta$lower <- rep(0, nrow(l_delta))
l_delta$median <- rep(NA, nrow(l_delta))
l_delta$upper <- rep(0, nrow(l_delta))

for(i in 1:nrow(l_delta)){
  l_delta[i, "median"] <- summary(get(paste0("l_fit_",unique_temp[which(l_delta[i,"temp"] == unique_temp_scaled)])))$coefficients["g",1]
}



EIP_plot <- EIP_percentile_plot_logistic(to_go_percentiles_plot, single_all_range_data, l_EIP)
delta_plot_ <- delta_plot_logistic(to_go_delta_percentiles_plot, delta_single_values, l_delta) +
  theme(legend.position = c(0.4,0.5))

legend <- cowplot::get_legend(delta_plot_)

x_axis <- ggdraw() + draw_label("Temperature / °C", color = "black", fontface = "plain", size = 20,
                                vjust = 0.01) + theme(plot.margin = margin(0, 0, 0, 0))
p <- plot_grid(EIP_plot + theme(axis.title.x = element_blank()), 
               delta_plot_ + theme(axis.title.x = element_blank(), legend.position = "none"))

p <- plot_grid(p, x_axis, ncol = 1, rel_heights = c(1, 0.05))

png(file = "figures/EIP_delta_plot.png", width = 1400, height = 450)
plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.2))
dev.off()

###########################################################################
################ mean parasite intensity functions ########################
###########################################################################
m_values <- generate_m_values(fit_extract, unique_temp_scaled, time_input, percentiles)

m_sensitivity_analysis <- m_values$m_sensitivity_analysis
m_EIP <- m_values$m_EIP

m_sensitivity_plot_ <- m_sensitivity_plot(m_sensitivity_analysis)
m_time_sensitivity_plot_ <- m_time_sensitivity_plot(m_EIP)

png(file = "figures/parasite_load_sensitivity_analysis.png", width = 1400, height = 500)
plot_grid(m_sensitivity_plot_ + theme(legend.key.size = unit(0.75,"cm"),
                                      legend.title = element_text(size = 18),
                                      legend.text = element_text(size = 20),
                                      legend.position = c(0.7, 0.25)), 
          m_time_sensitivity_plot_, ncol = 2)
dev.off()

################################################################################
################ prior and posterior distribution plots ########################
################################################################################
# plotting the posteriors
# calculating the mean number of parasites

parameters <- c("shape_oocyst", "m_rate_oocyst", "c_rate_oocyst",
                "shape_sporozoite", "rate_sporozoite", "mu_NB", "k_NB", 
                "m_delta", "c_delta", "a", "b", "beta_inf", "beta_temp", "sigma_error_delta", "sigma_error_survival") # "intercept", 
labels <- c("alpha[GO]", "m[beta]", "c[beta]", "alpha[OS]", "beta[OS]", "mu", "k", "m[delta]",
            "c[delta]", "a[s]", "b[s]", "beta[E]", "beta[C]", "sigma[delta]", "sigma[survival]") # "C[cox]",

posteriors <- parameter_extract_po(parameters, fit_new, labels)
posteriors$SD <- rep("2.5", nrow(posteriors))

posteriors_w <- parameter_extract_po(parameters, fit_wide, labels)
posteriors_w$SD <- rep("4.5", nrow(posteriors))

posteriors_all <- rbind(posteriors, posteriors_w)

a <- parameter_extract_pr(5, 25, 15.0, 2.5, "alpha[GO]")
b <- parameter_extract_pr( -5, 5, 0, 2.5, "m[beta]")
c <- parameter_extract_pr(0, 10, 3.0, 2.5, "c[beta]")
d <- parameter_extract_pr(5, 25, 15.0, 2.5, "alpha[OS]")
e <- parameter_extract_pr(0, 15, 1.875, 2.5, "beta[OS]")
f <- parameter_extract_pr(0, 7.5, 3.0, 2.5, "mu")
g <- parameter_extract_pr(0, 7.5, 0.1, 2.5, "k")
h <- parameter_extract_pr(-5, 5, 0, 2.5, "m[delta]")
i <- parameter_extract_pr(0, 10, 0.8, 2.5, "c[delta]")
j <- parameter_extract_pr(0, 0.25, 0.05, 2.5, "a[s]")
k <- parameter_extract_pr(0, 0.25, 0.05, 2.5, "b[s]")
l <- parameter_extract_pr(-4, 4, 0.5, 1.0, "beta[E]")
m <- parameter_extract_pr(-4, 4, 0.5, 1.0, "beta[C]")
#n <- parameter_extract_pr(-10, 5, -3, 2.5, "C[cox]")
o <- parameter_extract_pr(0, 1.5, 0, 0.25, "sigma[delta]")
p <- parameter_extract_pr(0, 1.5, 0, 0.25, "sigma[survival]")

priors <- rbind(a, b, c, d, e, f, g, h, i, j, k, l, m, o, p)

png(file = "figures/parameter_distribution_plot.png", width = 1250, height = 700)
ggplot() + geom_density(data= posteriors, aes(x = values, fill = as.factor(distribution), colour = as.factor(distribution)), 
                        stat = "density", alpha = 0.5) + 
  geom_area(data = priors, aes(x = x, y = y, 
                               fill = as.factor(distribution), colour = as.factor(distribution)), alpha = 0.5) + 
  facet_wrap(vars(label), scales = "free", labeller = label_parsed) + 
  scale_fill_manual(values = c("grey55", "grey25")) + scale_colour_manual(values = c("grey55", "grey25")) +
  theme(strip.text = element_text(face="bold", size=20)) +
  xlab("Parameter value") + ylab("Probability density") + blank_theme + 
  theme(legend.title = element_blank())
dev.off()

ggplot() + geom_density(data= posteriors_all, aes(x = values, fill = SD, colour = SD), 
                        stat = "density", alpha = 0.5) +
  facet_wrap(vars(label), scales = "free", labeller = label_parsed) + 
  theme(strip.text = element_text(face="bold", size=20)) +
  xlab("Parameter value") + ylab("Probability density") + blank_theme + 
  theme(legend.title = element_blank())

mean_new <- rstan::extract(fit_new, "shape_total_sporozoite[3]")$`shape_total_sporozoite[3]` /
  rstan::extract(fit_new, "rate_total_sporozoite[3]")$`rate_total_sporozoite[3]`

mean(mean_new)

mean_wide <- rstan::extract(fit_wide, "shape_total_sporozoite[3]")$`shape_total_sporozoite[3]` /
  rstan::extract(fit_wide, "rate_total_sporozoite[3]")$`rate_total_sporozoite[3]`

mean(mean_wide)

quantile(mean_new, c(0.05, 0.95))

quantile(mean_wide, c(0.05, 0.95))

quantile(rstan::extract(fit_new, "rate_oocyst[3]")$`rate_oocyst[3]`, c(0.05, 0.95))
quantile(rstan::extract(fit_new, "shape_oocyst")$`shape_oocyst`, c(0.05, 0.95))
quantile(rstan::extract(fit_new, "shape_sporozoite")$`shape_sporozoite`, c(0.05, 0.95))
quantile(rstan::extract(fit_new, "rate_sporozoite")$`rate_sporozoite`, c(0.05, 0.95))

#############################################################################################
################ extracting the EIP values at different temperatures ########################
#############################################################################################

EIP_table <- to_go_percentiles_plot[to_go_percentiles_plot$temp %in% unique_temp_scaled,]
for(i in 1:length(unique_temp)){
  EIP_table[which(EIP_table$temp == unique_temp_scaled[i]), "temp"] <- unique_temp[i]
}

EIP_table$value <- paste0(as.character(round(EIP_table$median, digits = 1))," (",
                          as.character(round(EIP_table$lower, digits = 1)),
                          " - ", as.character(round(EIP_table$upper, digits = 1)),")")
EIP_table <- EIP_table[,-c(3, 4, 5)]
EIP_table <- EIP_table %>% spread(EIP_X, value)

write.csv(EIP_table, "results/EIP_values_table.csv")

#########################################################################
################ extracting the parameter values ########################
#########################################################################

fit_summary <- summary(fit_new, pars = parameters)$summary
param_out <- as.data.frame(fit_summary)
param_out <- round(param_out, digits = 2)
write.csv(param_out, file = "results/parameter_values.csv")

#############################################################################
################ gamma distribution parameters at 27 ########################
#############################################################################
fit_27 <- readRDS(file = "results/mSOS_single_temperature_model_fit_new_temp_27")
fit_27_extract <- rstan::extract(fit_27)

prob_df <- generate_density_df(fit_27_extract)

png(file = "figures/individual_parasite_density_plot.png", width = 800, height = 300)
prob_plot(prob_df)
dev.off()

# extracting mean values
quantile(rstan::extract(fit_27, "shape_oocyst")$shape_oocyst/rstan::extract(fit_27, "rate_oocyst")$rate_oocyst, c(0.025, 0.5, 0.975))
quantile(rstan::extract(fit_27, "shape_sporozoite")$shape_sporozoite/rstan::extract(fit_27, "rate_sporozoite")$rate_sporozoite, c(0.025, 0.5, 0.975))
quantile(rstan::extract(fit_27, "shape_total_sporozoite")$shape_total_sporozoite/rstan::extract(fit_27, "rate_total_sporozoite")$rate_total_sporozoite, c(0.025, 0.5, 0.975))

# extracting the CDF values
subset(prob_df, density == "B (CDF)" & stage == "G to S" & lower >= 0.995)[1,"times"]
subset(prob_df, density == "B (CDF)" & stage == "G to S" & median >= 0.995)[1,"times"]
subset(prob_df, density == "B (CDF)" & stage == "G to S" & upper >= 0.995)[1,"times"]

#######################################################################
################ gamma distribution parameters ########################
#######################################################################

to_go_rate <- extract_rate(fit_extract, unique_temp)

to_go_rate_single_values <- extract_rate_single(indexes, unique_temp)

png(file = "figures/single_temp_model_param_plot.png", width = 800, height = 300)
single_temp_model_param_plot(to_go_rate_single_values, blank_theme)
dev.off()

