# R script 
# Author: Isaac J Stopard
# Version: 0.01 
# Last updated: 21/09/2020
# Notes: functions to put data in the correct format for plotting - note changes in packages may affect these functions

####################################################################
################ posterior predictive means ########################
####################################################################

# function that extracts the posterior predictive means from the MCMC output
prop_ppd_function <- function(fit_df, unique_temp, length_ppd_times, iterations, warmup, chains, length_unique_temp, Stan_data_name){
  prop_ppd <- array(NaN, c(length_ppd_times, ((iterations - warmup) * chains), length_unique_temp))
  for(i in 1:length_ppd_times){
    for(j in 1:length_unique_temp){
      prop_ppd[i,,j] <- fit_df[,paste0(Stan_data_name,"[",i,",",j,"]")]
    }
  }
  prop_ppd_df <- data.frame()
  labs <- c()
  for(i in 1:length_unique_temp){
    placeholder <- as.data.frame(prop_ppd[,,i])
    prop_ppd_df <- rbind(prop_ppd_df, placeholder)
    labs <- append(labs, rep(unique_temp[i], length_ppd_times))
  }
  day_post_inf <- rep(PPD_times,length_unique_temp)
  prop_ppd_df[,"temp"] <- labs
  prop_ppd_df[,"day_post_inf"] <- day_post_inf
  prop_ppd_df$median <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, median)
  prop_ppd_df$lower <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, quantile, probs = c(0.025))
  prop_ppd_df$upper <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, quantile, probs = c(0.975))
  prop_ppd_df$mean <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, mean)
  
  x <- ((iterations - warmup) * chains) 
  prop_ppd_df <- prop_ppd_df%>%gather("iteration", "value", 1:x)
  prop_ppd_df$percentile <- rep("iteration", nrow(prop_ppd_df))
  
  prop_quantile_ppd_df <- subset(prop_ppd_df, iteration == "V1")
  prop_quantile_ppd_df <- prop_quantile_ppd_df[,-c(7,8,9)]
  
  return(prop_quantile_ppd_df)
}


prop_ppd_function_single <- function(length_ppd_times, iterations, warmup, chains, fit, Stan_data_name){
  prop_ppd <- matrix(NaN, nrow = length_ppd_times, ncol = ((iterations - warmup) * chains))
  fit_df <- as.data.frame(fit)
  for(i in 1:length_ppd_times){
    prop_ppd[i,] <- fit_df[,paste0(Stan_data_name,"[",i,"]")]
  }
  prop_ppd_df <- cbind(as.data.frame(prop_ppd), as.data.frame(PPD_times))
  prop_ppd_df$median <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, median)
  prop_ppd_df$lower <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, quantile, probs = c(0.025))
  prop_ppd_df$upper <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, quantile, probs = c(0.975))
  prop_ppd_df$mean <- apply(prop_ppd_df[,1:((iterations - warmup)*chains)], 1, mean)
  
  x <- ((iterations - warmup) * chains) 
  prop_ppd_df <- prop_ppd_df%>%gather("iteration", "value", 1:x)
  
  prop_quantile_ppd_df <- subset(prop_ppd_df, iteration == "V1")
  prop_quantile_ppd_df <- prop_quantile_ppd_df[,-c(6,7)]
  
  colnames(prop_quantile_ppd_df) <- c("day_post_inf", "median", "lower", "upper", "mean")
  return(prop_quantile_ppd_df)
}

# generate the EIP for 27 degrees celsius
generate_EIP_27_data <- function(fit_27_extract, logistic_fit, label_){
  mSOS_EIP_df <- data.frame("times" = PPD_times, "lower" = rep(NaN, length(PPD_times)), "mSOS" = rep(NaN, length(PPD_times)), "upper" = rep(NaN, length(PPD_times)))
  for(i in 1:nrow(mSOS_EIP_df)){
    placeholder <- development_CDF_single(mSOS_EIP_df[i, "times"], fit_27_extract$shape_total_sporozoite, fit_27_extract$rate_total_sporozoite, fit_27_extract$mu_NB, fit_27_extract$k_NB)
    mSOS_EIP_df[i, "lower"] <- quantile(placeholder, 0.025)[[1]]
    mSOS_EIP_df[i, "mSOS"] <- quantile(placeholder, 0.5)[[1]]
    mSOS_EIP_df[i, "upper"] <- quantile(placeholder, 0.975)[[1]]
  }
  
  
  EIP_mSOS_lims <- cbind(mSOS_EIP_df, data.frame("logistic growth" = logistic_fit$median / summary(fitmodel27)$coefficients["g",1]))
  EIP_comparison <- gather(EIP_mSOS_lims[,c("times", "mSOS", "logistic.growth")], key = "model", value = "value", 2:3)
  EIP_mSOS_lims$label <- rep(label_, nrow(EIP_mSOS_lims))
  EIP_comparison$label <- rep(label_, nrow(EIP_comparison))
  return(list(EIP_mSOS_lims, EIP_comparison))
}

# plotting functions

oocyst_prevalence_plot <- function(oocyst_prop_ppd_df, oocyst_proportion_count_data, temp_, unique_temp_original, blank_theme){
  oocyst_prop_ppd_df$temp_label <- paste0(oocyst_prop_ppd_df$temp,"°C")
  oocyst_proportion_count_data$temp_label <- paste0(oocyst_proportion_count_data$temp,"°C")
  
  ggplot() +
    geom_ribbon(data = subset(oocyst_prop_ppd_df, temp == temp_ & day_post_inf <= 10), 
                aes(x = day_post_inf, ymin = lower, ymax = upper), fill = "grey55", alpha = 0.4) +
    geom_line(data = subset(oocyst_prop_ppd_df, temp == temp_ & day_post_inf <= 10), 
              aes(x = day_post_inf, y = median), colour = "black", linetype = 1, size = 0.75) +
    geom_pointrange(data = subset(oocyst_proportion_count_data, index_temp == which(unique_temp_original == temp_)), 
                    aes(x = day_post_inf, y = prevalence, ymin = lower, ymax = upper), colour = "black") +
    xlab("Days post blood feed") + ylab("Oocyst prevalence") + 
    scale_y_continuous(limits = c(0,0.95), breaks = seq(0,0.8,0.2), labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(limits = c(0,10), breaks = seq(0,10,2)) + blank_theme +
    facet_grid(. ~ temp_label) +  theme(strip.text = element_text(face="bold", size=20))
}

sporozoite_prevalence_plot <- function(sporozoite_prop_ppd_df, sporozoite_proportion_count_data, temp_, unique_temp_original, blank_theme){
  ggplot() +
    geom_ribbon(data = subset(sporozoite_prop_ppd_df, temp == temp_), 
                aes(x = day_post_inf, ymin = lower, ymax = upper), fill = "grey55", alpha = 0.4) +
    geom_line(data = subset(sporozoite_prop_ppd_df, temp == temp_), 
              aes(x = day_post_inf, y = median), colour = "black", linetype = 1, size = 0.75) +
    geom_pointrange(data = subset(sporozoite_proportion_count_data, index_temp == which(unique_temp_original == temp_)), 
                    aes(x = day_post_inf, y = prevalence, ymin = lower, ymax = upper), col = "black") +
    xlab("Days post blood feed") + ylab("Sporozoite prevalence / %") + 
    scale_y_continuous(limits = c(0,1.0), breaks = seq(0,1,0.2), labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(limits = c(0,30), breaks = seq(0,30,10)) + blank_theme + theme(axis.title=element_text(size = 12))
}

oocyst_intensity_plot <- function(oocyst_intensity_ppd_df, oocyst_intensity_data, temp_, blank_theme){
  oocyst_intensity_ppd_df$temp_label <- paste0(oocyst_intensity_ppd_df$temp, "°C")
  oocyst_intensity_data$temp_label <- paste0(oocyst_intensity_data$temp, "°C")
  ggplot() + 
    geom_boxplot(data = subset(oocyst_intensity_data, temp == temp_ & day_post_inf <= 10),
                 aes(x = day_post_inf, y = number_parasites, group = as.factor(day_post_inf))) +
    geom_ribbon(data = subset(oocyst_intensity_ppd_df, temp == temp_ & day_post_inf <= 10), 
                aes(x = day_post_inf, ymin = lower_m, ymax = upper_m), fill = "grey55", alpha = 0.4) +
    geom_line(data = subset(oocyst_intensity_ppd_df, temp == temp_ & day_post_inf <= 10), 
              aes(x = day_post_inf, y = median), colour = "black", linetype = 1, size = 0.75) +
    xlab("Days post blood feed") + ylab("Oocyst intensity") + 
    scale_y_continuous(limits = c(0,10), breaks = seq(0,10,2.5)) +
    scale_x_continuous(limits = c(0,11), breaks = seq(0,11,2)) + blank_theme +
    facet_grid(. ~ temp_label) +  theme(strip.text = element_text(face="bold", size=20))
}

# oocyst intensity plot where oocyst intensity is provided as a point range
oocyst_intensity_plot_PR <- function(oocyst_intensity_ppd_df, oocyst_intensity_summary, temp_, blank_theme){
  oocyst_intensity_summary$temp_label <- paste0(oocyst_intensity_summary$temp, "°C")
  oocyst_intensity_ppd_df$temp_label <- paste0(oocyst_intensity_ppd_df$temp, "°C")
  ggplot() + 
    geom_ribbon(data = subset(oocyst_intensity_ppd_df, temp == temp_ & day_post_inf <= 10), 
                aes(x = day_post_inf, ymin = lower_m, ymax = upper_m), fill = "grey55", alpha = 0.4) +
    geom_line(data = subset(oocyst_intensity_ppd_df, temp == temp_ & day_post_inf <= 10), 
              aes(x = day_post_inf, y = median), colour = "black", linetype = 1, size = 0.75) +
    geom_pointrange(data = subset(oocyst_intensity_summary, temp == temp_ & day_post_inf <= 10),
                    aes(x = day_post_inf, y = mean, ymin = lower, ymax = upper), col = "black") +
    xlab("Days post blood feed") + ylab("Oocyst intensity") + 
    scale_y_continuous(limits = c(0,10), breaks = seq(0,10,2.5)) +
    scale_x_continuous(limits = c(0,11), breaks = seq(0,11,2)) + blank_theme +
    facet_grid(. ~ temp_label) +  theme(strip.text = element_text(face="bold", size=20))
}


sporozoite_prevalence_plot_all <- function(sporozoite_prop_ppd_df, sporozoite_totals_all, blank_theme){
  sporozoite_prop_ppd_df$temp_label <- paste0(sporozoite_prop_ppd_df$temp, "°C")
  sporozoite_totals_all$temp_label <- paste0(sporozoite_totals_all$temp, "°C")
  ggplot() +
    geom_ribbon(data = sporozoite_prop_ppd_df, 
                aes(x = day_post_inf, ymin = lower, ymax = upper), fill = "grey55", alpha = 0.4) +
    geom_line(data = sporozoite_prop_ppd_df, 
              aes(x = day_post_inf, y = median), colour = "black", linetype = 1, size = 0.75) +
    geom_pointrange(data = sporozoite_totals_all, 
                    aes(x = day_post_inf, y = prevalence, ymin = lower, ymax = upper), col = "black") +
    facet_wrap(vars(temp_label), scales = "free") +
    xlab("Days post blood feed") + ylab("Sporozoite prevalence") + 
    scale_y_continuous(limits = c(0,0.95), breaks = seq(0,0.8,0.2), labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(limits = c(0,30), breaks = seq(0,30,10)) + blank_theme + 
    theme(strip.text = element_text(face="bold", size=20))
}

# function to plot the logistic fits
logistic_fit_plot_all <- function(sporozoite_totals_all, logistic_fit_df, unique_temp,
                                  day_cut_off, blank_theme){
  
  placeholder <- sporozoite_totals_all
  placeholder$included <- rep(NA, nrow(placeholder))
  single_temps <- unique_temp[-which(unique_temp == 33)]
  day_cut_off_temp <- day_cut_off[-which(unique_temp == 33)]
  for(i in 1:length(single_temps)){
    placeholder[which(placeholder$temp == single_temps[i] & placeholder$day_post_inf <= day_cut_off_temp[i]), "included"] <- "yes"
    placeholder[which(placeholder$temp == single_temps[i] & placeholder$day_post_inf > day_cut_off_temp[i]), "included"] <- "no"
  }
  
  placeholder$temp_label <- paste0(placeholder$temp, "°C")
  logistic_fit_df$temp_label <- paste0(logistic_fit_df$temp, "°C")
  
  ggplot() +
    geom_pointrange(data = subset(placeholder, temp!=33), 
                    aes(x = day_post_inf, y = prevalence, ymin = lower, ymax = upper, colour = as.factor(included))) +
    geom_line(data = logistic_fit_df, 
              aes(x = day_post_inf, y = median), colour = "black", linetype = 1, size = 0.75) +
    facet_wrap(vars(temp_label), scales = "free") +
    xlab("Days post blood feed") + ylab("Sporozoite prevalence") + 
    scale_y_continuous(limits = c(0,0.95), breaks = seq(0,0.8,0.2), labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(limits = c(0,30), breaks = seq(0,30,10)) + blank_theme + 
    scale_colour_manual(labels = c("No", "Yes"), values = c("grey55", "black")) +
    labs(colour = "Model fit \n to data") +
    theme(strip.text = element_text(face="bold", size=20),
          legend.key = element_rect(fill = "transparent", color = "transparent"), legend.text = element_text(size = 20),
          legend.background = element_rect(fill = "transparent", color = "transparent"))
}



plot_EIP_27 <- function(mSOS_EIP_df, EIP_comparison){
  ggplot() +
    geom_ribbon(data = subset(mSOS_EIP_df, times >= 5 & times <= 20), aes(x = times, ymin = lower, ymax = upper), fill = "grey55", alpha = 0.4) +
    geom_line(data = subset(EIP_comparison, times >= 5 & times <= 20), aes(x = times, y = value, colour = model), size = 0.85) +
    scale_colour_manual(labels = c("logistic growth", "mSOS"), values = c("black", "#56B4E9"), guide = guide_legend(reverse=TRUE)) +
    blank_theme + xlab("Days post blood feed") +
    ylab("EIP percentile") + 
    scale_x_continuous(limits = c(5,20)) + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    facet_grid(. ~ label) + 
    theme(strip.text = element_text(face="bold", size=20),
          legend.position = c(0.175, 0.875), legend.title = element_blank(),
          legend.key = element_rect(fill = "transparent", color = "transparent"), legend.text = element_text(size = 20),
          legend.background = element_rect(fill = "transparent", color = "transparent"))
}

# plotting the survival data as a Kaplan-Meier curve
# creates survival data frame that can be plot in ggplot
createSurvivalFrame <- function(f.survfit){
  f.frame <- NULL
  if(length(names(f.survfit$strata)) == 0){
    f.frame <- data.frame(time=f.survfit$time, n.risk=f.survfit$n.risk, n.event=f.survfit$n.event, n.censor = f.survfit
                          $n.censor, surv=f.survfit$surv, upper=f.survfit$upper, lower=f.survfit$lower)
    f.start <- data.frame(time=c(0, f.frame$time[1]), n.risk=c(f.survfit$n, f.survfit$n), n.event=c(0,0),
                          n.censor=c(0,0), surv=c(1,1), upper=c(1,1), lower=c(1,1))
    f.frame <- rbind(f.start, f.frame)
    f.frame$strata_infection <- rep(1, nrow(f.frame))
    rm(f.start)
  } else {
    f.strata <- NULL
    for(f.i in 1:length(f.survfit$strata)){
      f.strata <- c(f.strata, rep(names(f.survfit$strata)[f.i], f.survfit$strata[f.i]))
    }
    f.frame <- data.frame(time=f.survfit$time, n.risk=f.survfit$n.risk, n.event=f.survfit$n.event, n.censor = f.survfit
                          $n.censor, surv=f.survfit$surv, upper=f.survfit$upper, lower=f.survfit$lower, strata=factor(f.strata))
    rm(f.strata)
    for(f.i in 1:length(f.survfit$strata)){
      f.subset <- subset(f.frame, strata==names(f.survfit$strata)[f.i])
      f.start <- data.frame(time=c(0, f.subset$time[1]), n.risk=rep(f.survfit[f.i]$n, 2), n.event=c(0,0), 
                            n.censor=c(0,0), surv=c(1,1), upper=c(1,1), lower=c(1,1), strata=rep(names(f.survfit$strata)[f.i],
                                                                                                 2))
      f.frame <- rbind(f.start, f.frame)
      rm(f.start, f.subset)
    }
    f.frame <- f.frame[order(f.frame$strata, f.frame$time), ]
    rownames(f.frame) <- NULL
    f.frame$strata_infection <- rep(NaN, nrow(f.frame))
    f.frame[which(f.frame$strata == "initial_infection=1"),"strata_infection"] <- 1
    f.frame[which(f.frame$strata == "initial_infection=0"),"strata_infection"] <- 0
    f.frame <- f.frame[,-c(8)]
  }
  return(f.frame)
}

survival_plot <- function(surv_frame, St_ppd_df, blank_theme){
  surv_frame$temp_label <- paste0(surv_frame$temp, "°C")
  St_ppd_df$temp_label <- paste0(St_ppd_df$temp, "°C")
  ggplot(data=surv_frame) + 
    geom_ribbon(data = St_ppd_df, 
                aes(x = day_post_inf, ymin = lower, ymax = upper,
                    group = factor(infection_status)), alpha = 0.25, fill = "grey55") +
    #scale_fill_manual(labels = c("No", "Yes"), values = c("gray80","grey55")) +
    geom_line(data = St_ppd_df, 
              aes(x = day_post_inf, y = median, group = factor(infection_status), colour = factor(infection_status)), 
              size = 0.75) + #colour = "white"
    geom_step(aes(x=time, y=surv, group=factor(strata_infection)), 
              direction="hv", size = 1.5, colour="white") +
    geom_step(aes(x=time, y=surv, colour=factor(strata_infection), group=factor(strata_infection)),
              direction="hv", size = 1.2) +
    geom_step(aes(x=time, y=upper, colour=factor(strata_infection), group=factor(strata_infection))
              , direction="hv", linetype=2) + 
    geom_step(aes(x=time,y=lower, colour=factor(strata_infection), group=factor(strata_infection))
              , direction="hv", linetype=2) + 
    geom_point(data=subset(surv_frame, n.censor>=1), aes(x=time, y=surv, colour=factor(strata_infection), group=factor(strata_infection)),
               shape=3, size = 1.0, colour = "grey55") +
    scale_y_continuous(limits = c(0,1.0), breaks = seq(0,1,0.25), labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(limits = c(0,30), breaks = seq(0,30,10)) +
    ylab("Mosquito survival / percent alive") + xlab("Days post blood feed") +
    scale_colour_manual(labels = c("No", "Yes"), values = c("#E69F00", "#56B4E9")) +
    labs(colour = "Infectious \n blood fed", fill = "Infectious \n blood fed") + blank_theme + guides(colour = guide_legend(override.aes = list(size=7.5))) + 
    facet_wrap(vars(temp_label), scales = "free") + theme(strip.text = element_text(face="bold", size=20))
  
}


kaplan_meier_plot <- function(surv_frame, blank_theme){
  surv_frame$temp_label <- paste0(surv_frame$temp, "°C")
  ggplot(data=surv_frame) + 
    geom_step(aes(x=time, y=surv, group=factor(strata_infection)), 
              direction="hv", size = 1.5, colour="white") +
    geom_step(aes(x=time, y=surv, colour=factor(strata_infection), group=factor(strata_infection)),
              direction="hv", size = 1.2) +
    geom_step(aes(x=time, y=upper, colour=factor(strata_infection), group=factor(strata_infection))
              , direction="hv", linetype=2) + 
    geom_step(aes(x=time,y=lower, colour=factor(strata_infection), group=factor(strata_infection))
              , direction="hv", linetype=2) + 
    geom_point(data=subset(surv_frame, n.censor>=1), aes(x=time, y=surv, colour=factor(strata_infection), group=factor(strata_infection)),
               shape=3, size = 1.0, colour = "grey55") +
    scale_y_continuous(limits = c(0,1.0), breaks = seq(0,1,0.25), labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(limits = c(0,30), breaks = seq(0,30,10)) +
    ylab("Mosquito survival / percent alive") + xlab("Days post blood feed") +
    scale_colour_manual(labels = c("No", "Yes"), values = c("#E69F00", "#56B4E9")) +
    labs(colour = "Infectious \n blood fed", fill = "Infectious \n blood fed") + blank_theme + guides(colour = guide_legend(override.aes = list(size=7.5))) + 
    facet_wrap(vars(temp_label), scales = "free") + theme(strip.text = element_text(face="bold", size=20))
}


# plotting functions for single temperatures
sporozoite_prevalence_plot_single <- function(sporozoite_prop_ppd_df, sporozoite_proportion_count_data, temp_, blank_theme, label_){
  sporozoite_prop_ppd_df$label <- rep(label_, nrow(sporozoite_prop_ppd_df))
  sporozoite_proportion_count_data$label <- rep(label_, nrow(sporozoite_proportion_count_data))
  ggplot() +
    geom_ribbon(data = sporozoite_prop_ppd_df, 
                aes(x = day_post_inf, ymin = lower, ymax = upper), fill = "grey55", alpha = 0.4) +
    geom_line(data = sporozoite_prop_ppd_df, 
              aes(x = day_post_inf, y = median), colour = "black", linetype = 1, size = 0.75) +
    geom_pointrange(data = subset(sporozoite_proportion_count_data, temp == temp_), 
                    aes(x = day_post_inf, y = prevalence, ymin = lower, ymax = upper), col = "black") +
    xlab("Days post blood feed") + ylab("Sporozoite prevalence") + 
    scale_y_continuous(limits = c(0,0.95), breaks = seq(0,0.8,0.2), labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(limits = c(0,30), breaks = seq(0,30,10)) + blank_theme +
    facet_grid(. ~ label) +  theme(strip.text = element_text(face="bold", size=20))

}

sporozoite_prevalence_plot_logistic_single <- function(sporozoite_proportion_count_data, blank_theme, label_, temp_,
                                                       logistic_data){
  
  logistic_data$label <- rep(label_, nrow(logistic_data))
  sporozoite_data <- subset(sporozoite_proportion_count_data, temp == temp_)
  before <- subset(sporozoite_data, day_post_inf <= sporozoite_data[which(sporozoite_data$prevalence == max(sporozoite_data$prevalence)), "day_post_inf"])
  after <- subset(sporozoite_data, day_post_inf > sporozoite_data[which(sporozoite_data$prevalence == max(sporozoite_data$prevalence)), "day_post_inf"])
  before$included <- rep("yes", nrow(before))
  after$included <- rep("no", nrow(after))
  sporozoite_data <- rbind(before, after)
  sporozoite_data$label <- rep(label_, nrow(sporozoite_data))
  
  ggplot() +
    geom_line(data = logistic_data, 
              aes(x = day_post_inf, y = median), col = "black", size = 1.0) +
    geom_pointrange(data = sporozoite_data, 
                    aes(x = day_post_inf, y = prevalence, ymin = lower, ymax = upper, col = included)) +
    xlab("Days post blood feed") + ylab("Sporozoite prevalence") + 
    scale_y_continuous(limits = c(0,0.95), breaks = seq(0,0.8,0.2), labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(limits = c(0,30), breaks = seq(0,30,10)) +
    scale_colour_manual(labels = c("No", "Yes"), values = c("grey55", "black")) +
    blank_theme + facet_grid(. ~ label) +  labs(colour = "Data included \n in model fitting") +
    theme(strip.text = element_text(face="bold", size=20),
          legend.position = c(0.175, 0.875),
          legend.key = element_rect(fill = "transparent", color = "transparent"), legend.text = element_text(size = 20),
          legend.background = element_rect(fill = "transparent", color = "transparent"))
  
}

oocyst_prevalence_plot_single <- function(oocyst_prop_ppd_df, oocyst_proportion_count_data, temp_, blank_theme, label_){
  oocyst_prop_ppd_df$label <- rep(label_, nrow(oocyst_prop_ppd_df))
  oocyst_proportion_count_data$label <- rep(label_, nrow(oocyst_proportion_count_data))
  ggplot() +
    geom_ribbon(data = subset(oocyst_prop_ppd_df, day_post_inf <= 10), 
                aes(x = day_post_inf, ymin = lower, ymax = upper), fill = "grey55", alpha = 0.4) +
    geom_line(data = subset(oocyst_prop_ppd_df, day_post_inf <= 10), 
              aes(x = day_post_inf, y = median), colour = "black", linetype = 1, size = 0.75) +
    geom_pointrange(data = subset(oocyst_proportion_count_data, temp == temp_), 
                    aes(x = day_post_inf, y = prevalence, ymin = lower, ymax = upper), colour = "black") +
    xlab("Days post blood feed") + ylab("Oocyst prevalence") + 
    scale_y_continuous(limits = c(0,0.95), breaks = seq(0,0.8,0.2), labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(limits = c(0,10), breaks = seq(0,10,2)) + blank_theme +
    facet_grid(. ~ label) +  theme(strip.text = element_text(face="bold", size=20))
}

oocyst_intensity_plot_single <- function(oocyst_intensity_ppd_df, oocyst_intensity_data, temp_, blank_theme){
  ggplot() + 
    geom_boxplot(data = subset(oocyst_intensity_data, temp == temp_ & day_post_inf <= 10),
                 aes(x = day_post_inf, y = number_parasites, group = as.factor(day_post_inf))) +
    geom_ribbon(data = subset(oocyst_intensity_ppd_df, day_post_inf <= 10), 
                aes(x = day_post_inf, ymin = lower_m, ymax = upper_m), fill = "grey55", alpha = 0.4) +
    geom_line(data = subset(oocyst_intensity_ppd_df, day_post_inf <= 10), 
              aes(x = day_post_inf, y = median), colour = "black", linetype = 1, size = 0.75) +
    xlab("Days post blood feed") + ylab("Oocyst intensity") + 
    scale_y_continuous(limits = c(0,15), breaks = seq(0,15,5)) +
    scale_x_continuous(limits = c(0,11), breaks = seq(0,11,2)) + blank_theme
}


oocyst_intensity_plot_single_PR <- function(oocyst_intensity_ppd_df, oocyst_intensity_summary, temp_, blank_theme, label){
  oocyst_intensity_summary$label <- rep(label, nrow(oocyst_intensity_summary))
  oocyst_intensity_ppd_df$label <- rep(label, nrow(oocyst_intensity_ppd_df))
  ggplot() + 
    geom_ribbon(data = subset(oocyst_intensity_ppd_df, day_post_inf <= 10), 
                aes(x = day_post_inf, ymin = lower_m, ymax = upper_m), fill = "grey55", alpha = 0.4) +
    geom_line(data = subset(oocyst_intensity_ppd_df, day_post_inf <= 10), 
              aes(x = day_post_inf, y = median), colour = "black", linetype = 1, size = 0.75) +
    geom_pointrange(data = subset(oocyst_intensity_summary, temp == temp_ & day_post_inf <= 10),
                    aes(x = day_post_inf, y = mean, ymin = lower, ymax = upper), col = "black") +
    xlab("Days post blood feed") + ylab("Oocyst intensity") + 
    scale_y_continuous(limits = c(0,10), breaks = seq(0,10,2.5)) +
    scale_x_continuous(limits = c(0,11), breaks = seq(0,11,2)) + blank_theme +
    facet_grid(. ~ label) +  theme(strip.text = element_text(face="bold", size=20))
}


survival_plot_single <- function(f.frame, blank_theme, St_ppd_df, temp_, unique_temp, unique_temp_original){
  
  if("strata" %in% names(f.frame) == FALSE){
    ggplot(data=f.frame) + 
      geom_ribbon(data = St_ppd_df, 
                  aes(x = day_post_inf, ymin = lower, ymax = upper,
                      fill = factor(infection_status))) +
      scale_fill_manual(labels = c("No", "Yes"), values = c("gray80","grey55")) +
      geom_line(data = St_ppd_df,
                aes(x = day_post_inf, y = median, group = factor(infection_status)), 
                colour = "white", size = 0.4) +
      geom_step(aes(x=time, y=surv), direction="hv", size = 1.75, colour="white") +
      geom_step(aes(x=time, y=surv), direction="hv", size = 1.2, colour = "grey55") + 
      geom_step(aes(x=time, y=upper), direction="hv", linetype=2, colour = "grey55") + 
      geom_step(aes(x=time,y=lower), direction="hv", linetype=2, colour = "grey55") + 
      geom_point(data=subset(f.frame, n.censor>=1), aes(x=time, y=surv), shape=3, colour = "grey55") +
      ylab("Mosquito survival / percent alive") + xlab("Days post blood feed") +
      scale_y_continuous(limits = c(0,1.0), breaks = seq(0,1,0.2), labels = seq(0,1,0.2)*100) +
      scale_x_continuous(limits = c(0,30), breaks = seq(0,30,5)) + 
      labs(fill = "Infectious \n blood fed") + blank_theme + guides(colour = guide_legend(override.aes = list(size=7.5)))
  } else{
    ggplot(data=f.frame) + 
      geom_ribbon(data = St_ppd_df, 
                  aes(x = day_post_inf, ymin = lower, ymax = upper,
                      fill = factor(infection_status))) +
      scale_fill_manual(labels = c("No", "Yes"), values = c("gray80","grey55")) +
      geom_line(data = St_ppd_df, 
                aes(x = day_post_inf, y = median, group = factor(infection_status)), 
                colour = "white", size = 0.4) +
      geom_step(aes(x=time, y=surv, colour=factor(strata_infection), group=factor(strata_infection)), 
                direction="hv", size = 1.75, colour="white") +
      geom_step(aes(x=time, y=surv, colour=factor(strata_infection), group=factor(strata_infection)),
                direction="hv", size = 1.2) +
      geom_step(aes(x=time, y=upper, colour=factor(strata_infection), group=factor(strata_infection))
                , direction="hv", linetype=2) + 
      geom_step(aes(x=time,y=lower, colour=factor(strata_infection), group=factor(strata_infection))
                , direction="hv", linetype=2) + 
      geom_point(data=subset(f.frame, n.censor>=1), aes(x=time, y=surv, colour=factor(strata_infection), group=factor(strata_infection)),
                 shape=3, size = 1.0, colour = "grey55") +
      scale_y_continuous(limits = c(0,1.0), breaks = seq(0,1,0.2), labels = seq(0,1,0.2)*100) +
      scale_x_continuous(limits = c(0,30), breaks = seq(0,30,5)) +
      ylab("Mosquito survival / percent alive") + xlab("Days post blood feed") +
      scale_colour_manual(labels = c("No", "Yes"), values = c("gray80","grey55")) +
      labs(colour = "Infectious \n blood fed") + blank_theme + guides(colour = guide_legend(override.aes = list(size=7.5)))
  }
}


##################################################################
################ EIP estimation functions ########################
##################################################################

# function that calculates the first order CDF of the all temperature model
development_CDF <- function(times, temps, shape_oocyst, m_rate_oocyst, c_rate_oocyst, 
                            shape_sporozoite, rate_sporozoite, mu, k){
  rate_oocyst <- (m_rate_oocyst * temps) + c_rate_oocyst
  mean_sporozoite <- (rate_oocyst * shape_sporozoite + shape_oocyst * rate_sporozoite)/(rate_oocyst * rate_sporozoite)
  variance_sporozoite <- (rate_oocyst^2 * shape_sporozoite + rate_sporozoite^2 * shape_oocyst) / (rate_oocyst^2 * rate_sporozoite^2)
  shape_total_sporozoite <- mean_sporozoite^2 / variance_sporozoite
  rate_total_sporozoite <- mean_sporozoite / variance_sporozoite
  regularised_gamma <- Rgamma(shape_total_sporozoite, (times * rate_total_sporozoite))
  parasite_CDF <- (1 - (k^k) * ((k + mu * regularised_gamma)^-k)) / (1 - (k / (k + mu))^k)
  return(parasite_CDF)
}

# function that returns the time at which the CDF is closest
EIP_time_function <- function(i, temp, p, parameter_df, times){
  a <- development_CDF(times, temp, parameter_df[i,"shape_oocyst"], parameter_df[i,"m_rate_oocyst"], 
                       parameter_df[i,"c_rate_oocyst"], parameter_df[i,"shape_sporozoite"], 
                       parameter_df[i,"rate_sporozoite"], parameter_df[i,"mu"], 
                       parameter_df[i,"k"])
  EIP <- times[which(abs(a-p) == min(abs(a-p)))]
  return(EIP)
}

generate_EIP_values <- function(parameter_df, time_input, temp_input, unique_temp_scaled, iteration_index, percentiles){
  to_go <- expand.grid(iteration_index, temp_input, percentiles)
  colnames(to_go) <- c("iteration_index","temps","EIP_percentile")
  EIP_times <- apply(to_go, 1, function(to_go) EIP_time_function(to_go[1],to_go[2], to_go[3], parameter_df = parameter_df, times = time_input)) # takes a long time to run
  to_go <- cbind(to_go, as.data.frame(EIP_times))
  colnames(to_go) <- c("iteration_index","temps","EIP_percentile","EIP_times")
  saveRDS(to_go, file = "results/hierarchical_mSOS_EIP_values")
}

generate_EIP_quantiles <- function(to_go, temps, percentiles, quantiles){
  to_go_percentiles <- expand.grid(temps, percentiles, quantiles)
  colnames(to_go_percentiles) <- c("temp","EIP_X", "quantile")
  to_go_percentiles$value <- rep(NaN, nrow(to_go_percentiles))
  for(i in 1:nrow(to_go_percentiles)){
    a <- subset(to_go, temps == to_go_percentiles[i,"temp"] &
                  EIP_percentile == to_go_percentiles[i,"EIP_X"])
  
    to_go_percentiles[i,"value"] <- quantile(a$EIP_times, probs = to_go_percentiles[i,"quantile"])[[1]]
  }

  to_go_percentiles_plot <- to_go_percentiles%>%spread(quantile, value)
  colnames(to_go_percentiles_plot) <- c("temp", "EIP_X","lower","median","upper")
  return(to_go_percentiles_plot)
}

# when run for each temperature independentily
development_CDF_single <- function(times, shape_sporozoite, rate_sporozoite, mu, k){
  regularised_gamma <- Rgamma(shape_sporozoite, (times * rate_sporozoite))
  parasite_CDF <- (1 - (k^k) * ((k + mu * regularised_gamma)^-k)) / (1 - (k / (k + mu))^k)
  # parasite_CDF <- 1 - (k / (k + mu * regularised_gamma))^k
  # parasite_CDF <- 1 + (-1 + regularised_gamma) * (k / (k + mu * regularised_gamma))^k
  return(parasite_CDF)
}

EIP_time_function_single <- function(i, p, parameter_df, times){
  a <- development_CDF_single(times, parameter_df[i,"shape_sporozoite"], 
                              parameter_df[i,"rate_sporozoite"], parameter_df[i,"mu"], 
                              parameter_df[i,"k"])
  EIP <- times[which(abs(a-p) == min(abs(a-p)))]
  return(EIP)
}

generate_EIP_values_single <- function(unique_temp_original, iteration_index_single, times, percentiles){
  to_go_single <- expand.grid(iteration_index_single, percentiles)
  colnames(to_go_single) <- c("iteration_index_single", "percentile")
  EIP_times_single_all <- c()
  temps_single_all <- c()
  temp_indexes <- c(1,2,3,4,5,7)
  for(i in 1:length(temp_indexes)){
    index <- temp_indexes[i]
    if(unique_temp_original[index] == 27){
      fit_single <- readRDS(file = "results/mSOS_single_temperature_model_fit_new_temp_27")
      fit_extract_single <- rstan::extract(fit_single)
      parameter_df_single <- data.frame(iteration = iteration_index_single,
                                        shape_sporozoite = fit_extract_single$shape_total_sporozoite,
                                        rate_sporozoite = fit_extract_single$rate_total_sporozoite,
                                        mu = fit_extract_single$mu_NB,
                                        k = fit_extract_single$k_NB)
      
      EIP_times_single <- rep(NaN, nrow(to_go_single))
      temps_single <- rep(round(unique_temp[index], digits = 9), nrow(to_go_single))
      for(j in 1:nrow(to_go_single)){
        EIP_times_single[j] <- EIP_time_function_single(to_go_single[j,1], to_go_single[j,2],
                                                        parameter_df = parameter_df_single, times = times)
      }
      
      EIP_times_single_all <- append(EIP_times_single_all, EIP_times_single)
      temps_single_all <- append(temps_single_all, temps_single)
    } else{
      fit_single <- readRDS(file = paste0("results/mSOS_single_temperature_model_fit_new_temp_",unique_temp_original[index]))
      fit_extract_single <- rstan::extract(fit_single)
      parameter_df_single <- data.frame(iteration = iteration_index_single,
                                        shape_sporozoite = fit_extract_single$shape_sporozoite,
                                        rate_sporozoite = fit_extract_single$rate_sporozoite,
                                        mu = fit_extract_single$mu_NB,
                                        k = fit_extract_single$k_NB)
      
      EIP_times_single <- rep(NaN, nrow(to_go_single))
      temps_single <- rep(unique_temp[index], nrow(to_go_single))
      
      for(j in 1:nrow(to_go_single)){
        EIP_times_single[j] <- EIP_time_function_single(to_go_single[j,1], to_go_single[j,2],
                                                        parameter_df = parameter_df_single, times = times)
      }
      
      EIP_times_single_all <- append(EIP_times_single_all, EIP_times_single)
      temps_single_all <- append(temps_single_all, temps_single)
    }
  } 
  
  b <- rbind(to_go_single, to_go_single, to_go_single, to_go_single, to_go_single, to_go_single)
  c <- cbind(as.data.frame(EIP_times_single_all), as.data.frame(temps_single_all))
  to_go_single_all <- cbind(b, c)
  saveRDS(to_go_single_all, file = "results/mSOS_EIP_values_single")
  return(to_go_single_all)
}

EIP_percentile_plot <- function(to_go_percentiles_plot, single_all_range_data){ 
  single_all_range_data$label <- rep("A", nrow(single_all_range_data))
  to_go_percentiles_plot$label <- rep("A", nrow(to_go_percentiles_plot))
  
  ggplot() +
  geom_ribbon(data = to_go_percentiles_plot, aes(x = temp, ymin = lower, ymax = upper, group = factor(EIP_X)), fill = "grey55", alpha = 0.4) +
  geom_line(data = to_go_percentiles_plot, aes(x = temp, y = median, colour = factor(EIP_X)), size = 0.75) +
  geom_pointrange(data = single_all_range_data, aes(x = temp, y = median, ymin = lower, ymax = upper, colour = factor(EIP_X))) +
  scale_color_manual(labels = c(parse(text="EIP[10]"),parse(text="EIP[50]"),parse(text="EIP[90]")), values = c("#E69F00", "#56B4E9", "#CC79A7")) + 
  #scale_fill_manual(labels = c(parse(text="EIP[10]"),parse(text="EIP[50]"),parse(text="EIP[90]")), values = c("grey35", "grey55", "grey75")) +
  #scale_shape_manual(labels = c(parse(text="EIP[10]"),parse(text="EIP[50]"),parse(text="EIP[90]")), values = c(16,15,17)) +
  theme(legend.position = c(0.75, 0.75)) + blank_theme + theme(legend.text = element_text(size = 15)) + 
  scale_x_continuous(limits = c(unique_temp_scaled[1], unique_temp_scaled[7]), breaks = seq(unique_temp_scaled[1],unique_temp_scaled[7], 
                                                                                            (unique_temp_scaled[7]-unique_temp_scaled[6])*3), 
                     labels = seq(unique_temp[1], unique_temp[7], 3)) + xlab("Temperature / °C") +
  scale_y_continuous(limits = c(4,30), breaks = seq(5,30,5)) + ylab("EIP / days") +
  labs(colour = "", shape = "", fill = "") + facet_wrap(.~label, scales = "free") +
  theme(strip.text = element_text(face="bold", size=20)) # + guides(colour = FALSE) 
}

# plot that includes the logistic model EIP points
EIP_percentile_plot_logistic <- function(to_go_percentiles_plot, single_all_range_data, l_EIP){ 
  
  l_EIP$model <- rep("logistic growth", nrow(l_EIP))
  single_all_range_data$model <- rep("mSOS", nrow(single_all_range_data))
  single_all_range_data <- rbind(l_EIP, single_all_range_data)
  to_go_percentiles_plot$model <- rep("all temp mSOS", nrow(to_go_percentiles_plot))
  single_all_range_data$label <- rep("A", nrow(single_all_range_data))
  to_go_percentiles_plot$label <- rep("A", nrow(to_go_percentiles_plot))
  single_for_legend <- rbind(single_all_range_data, to_go_percentiles_plot)
 
  ggplot() +
    geom_ribbon(data = to_go_percentiles_plot, aes(x = temp, ymin = lower, ymax = upper, group = factor(EIP_X)), fill = "grey55", alpha = 0.4) +
    geom_line(data = to_go_percentiles_plot, aes(x = temp, y = median, colour = factor(EIP_X)), size = 1) +
    geom_linerange(data = subset(single_all_range_data, model == "mSOS"), aes(x = temp, ymin = lower, ymax = upper, 
                                                                              colour = factor(EIP_X)), size = 1) +
    geom_point(data = single_all_range_data, aes(x = temp, y = median, fill = factor(EIP_X), shape = model), size = 4, alpha = 0.75) +
    scale_color_manual(labels = c(parse(text="EIP[10]"),parse(text="EIP[50]"),parse(text="EIP[90]")), values = c("#E69F00", "#56B4E9", "#CC79A7")) + 
    scale_fill_manual(labels = c(parse(text="EIP[10]"),parse(text="EIP[50]"),parse(text="EIP[90]")), values = c("#E69F00", "#56B4E9", "#CC79A7")) +
    scale_shape_manual(breaks = c("logistic growth", "mSOS", "all temp mSOS"), labels = c("logistic growth", "mSOS (single temp)", "mSOS (all temp)"), 
                       values = c(24, 21, 95)) +
    theme(legend.position = c(0.9, 0.9)) + blank_theme + 
    scale_x_continuous(limits = c(unique_temp_scaled[1], unique_temp_scaled[7]), breaks = seq(unique_temp_scaled[1],unique_temp_scaled[7], 
                                                                                              (unique_temp_scaled[7]-unique_temp_scaled[6])*3), 
                       labels = seq(unique_temp[1], unique_temp[7], 3)) + xlab("Temperature / °C") +
    scale_y_continuous(limits = c(4,27.5), breaks = seq(5,25,5)) + ylab("EIP / days") +
    labs(colour = "", shape = "", fill = "") + facet_wrap(.~label, scales = "free") +
    theme(strip.text = element_text(face="bold", size=20), legend.key = element_rect(fill = "transparent", color = "transparent"),
          legend.background = element_rect(fill = "transparent", color = "transparent"),
          legend.box = "horizontal", legend.text = element_text(size = 17)) +
    geom_blank(data = single_for_legend, aes(x = temp, y = median, shape = factor(model))) +
    guides(colour = FALSE, shape = FALSE,
           fill = guide_legend(override.aes = list(shape = c(15, 15, 15), colour = c("#E69F00", "#56B4E9", "#CC79A7"))))
}


#########################################################
################ delta functions ########################
#########################################################

extract_delta <- function(fit_extract, iteration_index, temps, quantiles){
  parameter_delta <- data.frame(iteration = seq(1:(iterations - warmup) * chains),
                                m_delta = fit_extract$m_delta,
                                c_delta = fit_extract$c_delta)

  to_go_delta <- expand.grid(iteration_index, temps)
  colnames(to_go_delta) <- c("iteration", "temp")
  to_go_delta$delta <- rep(NaN, nrow(to_go_delta))

  # generating the delta values for each iteration
  for(i in 1:nrow(to_go_delta)){
    to_go_delta[i,"delta"] = 1 / (1 + exp(
      -((parameter_delta[to_go_delta[i,"iteration"],"m_delta"] * to_go_delta[i,"temp"]) +
          parameter_delta[to_go_delta[i, "iteration"], "c_delta"])))
  }

  to_go_delta_percentiles <- expand.grid(temps, quantiles)
  colnames(to_go_delta_percentiles) <- c("temp", "percentile")
  to_go_delta_percentiles$value <- rep(NaN, nrow(to_go_delta_percentiles))
  for(i in 1:nrow(to_go_delta_percentiles)){
    a <- subset(to_go_delta, temp == to_go_delta_percentiles[i,"temp"])
    to_go_delta_percentiles[i,"value"] <- quantile(a$delta, to_go_delta_percentiles[i,"percentile"])[[1]]
  }

  to_go_delta_percentiles_plot <- to_go_delta_percentiles%>%spread(percentile, value)
  colnames(to_go_delta_percentiles_plot) <- c("temp", "lower", "median", "upper")
  return(to_go_delta_percentiles_plot)
}


extract_delta_single <- function(temp_indexes, iteration_index_single, unique_temp, unique_temp_original){
  
  delta_single_values <- as.data.frame(unique_temp[temp_indexes])
  colnames(delta_single_values) <- c("temp")
  delta_single_values$lower <- rep(NaN, nrow(delta_single_values))
  delta_single_values$median <- rep(NaN, nrow(delta_single_values))
  delta_single_values$upper <- rep(NaN, nrow(delta_single_values))

  for(i in 1:length(temp_indexes)){
    index <- temp_indexes[i]
    fit_single_delta <- readRDS(file = paste0("results/mSOS_single_temperature_model_fit_new_temp_",unique_temp_original[index]))
    fit_extract_single_delta <- rstan::extract(fit_single_delta)
    parameter_delta_single <- data.frame(delta = fit_extract_single_delta$delta)
    delta_single_values[i,"lower"] <- quantile(parameter_delta_single[,"delta"], 0.025)
    delta_single_values[i,"median"] <- quantile(parameter_delta_single[,"delta"], 0.5)
    delta_single_values[i,"upper"] <- quantile(parameter_delta_single[,"delta"], 0.975)
  }
  return(delta_single_values)
}


delta_plot <- function(to_go_delta_percentiles_plot, delta_single_values){
  to_go_delta_percentiles_plot$label <- rep("B", nrow(to_go_delta_percentiles_plot))
  delta_single_values$label <- rep("B", nrow(delta_single_values))
  
  ggplot() +
  #geom_line(data = to_go_delta, aes(x = temp, y = delta, group = factor(iteration)), colour = "grey", size = 0.05, alpha = 0.15) +
  geom_ribbon(data = to_go_delta_percentiles_plot, aes(x = temp, ymin = lower, ymax = upper), fill = "grey55", alpha = 0.4) +
  geom_line(data = to_go_delta_percentiles_plot, aes(x = temp, y = median)) +
  geom_pointrange(data=delta_single_values, aes(x = temp, y = median, ymin = lower, ymax = upper), colour = "black") +
  blank_theme +
  #scale_colour_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"), labels = c("5%","50%","95%")) +
  scale_x_continuous(limits = c(unique_temp_scaled[1], unique_temp_scaled[7]), breaks = seq(unique_temp_scaled[1],
                                                                                            unique_temp_scaled[7], (unique_temp_scaled[7]-unique_temp_scaled[6])*3), 
                     labels = seq(unique_temp[1], unique_temp[7], 3)) + xlab("Temperature / °C") + 
  ylab(expression(paste("Human-to-mosquito transmission probability (",delta,")"))) +
  scale_y_continuous(limits = c(0.3, 1.0), breaks = seq(0.3,1.0,0.1)) +
  facet_wrap(.~label, scales = "free") +
  theme(strip.text = element_text(face="bold", size=20), axis.title.y = element_text(size = 18))
}

delta_plot_logistic <- function(to_go_delta_percentiles_plot, delta_single_values, l_delta){
  l_delta$label <- rep("B", nrow(l_delta))
  l_delta$model <- rep("logistic growth", nrow(l_delta))
  
  to_go_delta_percentiles_plot$label <- rep("B", nrow(to_go_delta_percentiles_plot))
  to_go_delta_percentiles_plot$model <- rep("mSOS (all temp)", nrow(to_go_delta_percentiles_plot))
  
  delta_single_values$label <- rep("B", nrow(delta_single_values))
  delta_single_values$model <- rep("mSOS (single temp)", nrow(delta_single_values))
  
  delta_for_legend <- rbind(l_delta, to_go_delta_percentiles_plot, delta_single_values)
  
  ggplot() +
    geom_ribbon(data = to_go_delta_percentiles_plot, aes(x = temp, ymin = lower, ymax = upper), fill = "grey55", alpha = 0.4) +
    geom_line(data = to_go_delta_percentiles_plot, aes(x = temp, y = median)) +
    geom_pointrange(data=delta_single_values, aes(x = temp, y = median, ymin = lower, ymax = upper), colour = "black", size = 0.8, alpha = 0.9) +
    geom_point(data = l_delta, aes(x = temp, y = median, shape = factor(model)), colour = "black", size = 4.5, alpha = 0.9) +
    geom_blank(data = delta_for_legend, aes(x = temp, y = median, shape = factor(model))) +
    blank_theme +
    scale_shape_manual(breaks = c("logistic growth", "mSOS (single temp)", "mSOS (all temp)"), values = c(17, 95, 16)) +
    scale_x_continuous(limits = c(unique_temp_scaled[1], unique_temp_scaled[7]), breaks = seq(unique_temp_scaled[1],
                                                                                              unique_temp_scaled[7], (unique_temp_scaled[7]-unique_temp_scaled[6])*3), 
                       labels = seq(unique_temp[1], unique_temp[7], 3)) + xlab("Temperature / °C") + 
    ylab(expression(paste("Human-to-mosquito transmission probability (",delta,")"))) + 
    labs(shape = "") + scale_y_continuous(limits = c(0.3, 1.0), breaks = seq(0.3,1.0,0.1)) +
    facet_wrap(.~label, scales = "free") +
    theme(strip.text = element_text(face="bold", size=20), axis.title.y = element_text(size = 18)) +
    guides(shape = guide_legend(override.aes = list(size = c(4.5, 4.5, 10), shape = c(2, 1, 95)))) +
    theme(strip.text = element_text(face="bold", size=20), legend.key = element_rect(fill = "transparent", color = "transparent"),
          legend.background = element_rect(fill = "transparent", color = "transparent"),
          legend.box = "horizontal", legend.text = element_text(size = 17), legend.position = c(0.8, 0.9))
}

###########################################################################
################ mean parasite intensity functions ########################
###########################################################################

# m sensitivity analysis
# survival functions
cox <- function(beta_inf, beta_temp, infec, temp){
  result <- (beta_inf * infec) + (beta_temp * temp);
  return(result)
}

baseline_hazard <- function(t, a, b){
  result <- log(a) + (b * t);
  return(result)
}

S_t <- function(t, a, b, cox_){
  result <- (a/b)*(1-exp(b*t))*exp(cox_);
  return(result)
}

m_development_CDF <- function(times, shape_total_sporozoite, rate_total_sporozoite, mu, k){
  regularised_gamma <- Rgamma(shape_total_sporozoite, (times * rate_total_sporozoite))
  parasite_CDF <- (1 - (k^k) * ((k + mu * regularised_gamma)^-k)) / (1 - (k / (k + mu))^k)
  return(parasite_CDF)
}


generate_m_values <- function(fit_extract, unique_temp, times, percentiles){
  m_shape_oocyst <- mean(fit_extract$shape_oocyst)
  m_m_rate_oocyst <- mean(fit_extract$m_rate_oocyst)
  m_c_rate_oocyst <- mean(fit_extract$c_rate_oocyst)
  m_shape_sporozoite <- mean(fit_extract$shape_sporozoite)
  m_rate_sporozoite <- mean(fit_extract$rate_sporozoite)
  m_m_delta <- mean(fit_extract$m_delta)
  m_c_delta <- mean(fit_extract$c_delta)
  m_k_NB <- mean(fit_extract$k_NB)
  m_a <- mean(fit_extract$a)
  m_b <- mean(fit_extract$b)
  m_beta_inf <- mean(fit_extract$beta_inf)
  m_beta_temp <- mean(fit_extract$beta_temp)
  m_times <- seq(0,30,0.1)

  S_t_u <- exp(S_t(m_times, m_a, m_b, cox(m_beta_inf, m_beta_temp, 0, unique_temp[3])))
  S_t_i <- exp(S_t(m_times, m_a, m_b, cox(m_beta_inf, m_beta_temp, 1, unique_temp[3])))
  RR <- S_t_i / S_t_u
  m_delta <- (m_m_delta * unique_temp[3]) + m_c_delta

  m_rate_oocyst <- (m_m_rate_oocyst * unique_temp[3]) + m_c_rate_oocyst
  m_mean_sporozoite <- (m_rate_oocyst * m_shape_sporozoite + m_shape_oocyst * m_rate_sporozoite)/(m_rate_oocyst * m_rate_sporozoite)
  m_variance_sporozoite <- (m_rate_oocyst^2 * m_shape_sporozoite + m_rate_sporozoite^2 * m_shape_oocyst) / (m_rate_oocyst^2 * m_rate_sporozoite^2);
  m_shape_total_sporozoite <- m_mean_sporozoite^2 / m_variance_sporozoite
  m_rate_total_sporozoite <- m_mean_sporozoite / m_variance_sporozoite

  m <- seq(1,25,2)
  
  m_sensitivity_analysis <- expand.grid(m_times, m)
  colnames(m_sensitivity_analysis) <- c("time", "m")
  
  prevalence <- c()
  for(i in 1:length(m)){
    m_CDF <- m_development_CDF(m_times, m_shape_total_sporozoite, m_rate_total_sporozoite, m[i], m_k_NB)
    m_prevalence <- (m_CDF * RR) * m_delta
    prevalence <- append(prevalence, m_prevalence)
  }

  m_sensitivity_analysis <- cbind(m_sensitivity_analysis, prevalence)
  m_sensitivity_analysis$mean_parasite_intensity <- m_sensitivity_analysis$m / (1 - (m_k_NB/(m_sensitivity_analysis$m + m_k_NB))^m_k_NB)
  m_sensitivity_analysis[which(m_sensitivity_analysis$prevalence <0), "prevalence"] <- 0
  
  m_EIP_vals <- c(seq(1,10, 0.5),seq(11,25,1))
  
  m_EIP <- expand.grid(m_EIP_vals, percentiles)
  colnames(m_EIP) <- c("m", "percentile")
  m_EIP$time <- rep(NaN, nrow(m_EIP))
  for(i in 1:nrow(m_EIP)){
    m_CDF <- development_CDF(times, unique_temp[3], m_shape_oocyst, m_m_rate_oocyst, m_c_rate_oocyst, 
                             m_shape_sporozoite, m_rate_sporozoite, m_EIP[i,"m"], m_k_NB)
    m_EIP[i,"time"] <- times[which(abs(m_CDF - m_EIP[i,"percentile"]) == min(abs(m_CDF - m_EIP[i,"percentile"])))]
  }
  
  m_EIP$mean_parasite_intensity <- m_EIP$m / (1 - (m_k_NB/(m_EIP$m + m_k_NB))^m_k_NB)
  
  out <- list("m_sensitivity_analysis" = m_sensitivity_analysis, "m_EIP" = m_EIP)
  return(out)
}

# "#E69F00", "#56B4E9", "#CC79A7"

m_sensitivity_plot <- function(m_sensitivity_analysis){
  m_sensitivity_analysis$label <- rep("A", nrow(m_sensitivity_analysis))
  ggplot() +
    geom_line(data = m_sensitivity_analysis, aes(x = time, y = prevalence, group = mean_parasite_intensity, 
                                               colour = mean_parasite_intensity), size = 0.8) +
    xlab("Days post blood feed") + ylab("Sporozoite prevalence") + 
    scale_y_continuous(limits = c(0,0.675), breaks = seq(0,0.6,0.15), labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(limits = c(0,30), breaks = seq(0,30,5)) + 
    labs(colour = c("Infected mosquito\nmean parasite load")) + scale_color_gradient(low = "#56B4E9", high = "#CC79A7") +
    blank_theme + theme(legend.text=element_text(size=14), legend.title = element_text(size = 14), legend.position = c(0.85, 0.8)) +
    facet_wrap(.~label) + theme(strip.text = element_text(face="bold", size=20))
}

m_time_sensitivity_plot <- function(m_EIP){
  m_EIP$label_ <- rep("B", nrow(m_EIP))
  labels_df <- data.frame(x = c(24, 24, 24), y = c(7.25, 8.5, 10),
                          text = c("EIP[10]", "EIP[50]", "EIP[90]"), label_ = c("B", "B", "B"))
  ggplot() +
    geom_line(data = m_EIP, aes(x = mean_parasite_intensity, y = time, group = factor(percentile)),size = 1.5) +
    #scale_color_manual(labels = c(parse(text="EIP[10]"),parse(text="EIP[50]"),parse(text="EIP[90]")), 
    #                   values = c("grey35", "grey55", "grey75")) + 
    blank_theme + xlab("Infected mosquito mean parasite load") + ylab("EIP / days") +
    labs(colour = "") +  theme(legend.position = "none") +
    geom_text(data = labels_df, aes(x = x, y = y, label = text), parse = TRUE, size = 7) +
    facet_wrap(.~label_) + theme(strip.text = element_text(face="bold", size=20))
}

################################################################################
################ prior and posterior distribution plots ########################
################################################################################

parameter_extract_po <- function(parameters_, fit, labels_){
df_out <- data.frame("values" = NULL, "distribution" = NULL, "label" = NULL)
  for(i in 1:length(parameters_)){
    parameter <- parameters_[i]
    label <- labels_[i]
    df<- as.data.frame(rstan::extract(fit, parameter))
    df$distribution <- rep("posterior", nrow(df))
    df$label <- rep(label, nrow(df))
    colnames(df) <- c("values", "distribution", "label")
    df_out <- rbind(df_out, df)
  }
  return(df_out)
}

parameter_extract_pr <- function(lower, upper, mean, sd, label){
  parameter_range <- seq(lower, upper, 0.01)
  prior <- dnorm(parameter_range, mean, sd)
  prior_df <- cbind(as.data.frame(parameter_range), as.data.frame(prior))
  prior_df$distribution <- rep("prior", nrow(prior_df))
  colnames(prior_df) <- c("x", "y", "distribution")
  prior_df$label <- rep(label, nrow(prior_df))
  return(prior_df)
}

###############################################################
################ oocyst rate parameter ########################
###############################################################

extract_rate <- function(fit_extract, unique_temp){
  values <- as.data.frame(unique_temp)
  colnames(values) <- c("temp")
  
  values$lower_rate <- rep(NaN, nrow(values))
  values$median_rate <- rep(NaN, nrow(values))
  values$upper_rate <- rep(NaN, nrow(values))
  
  values$lower_shape <- rep(NaN, nrow(values))
  values$median_shape <- rep(NaN, nrow(values))
  values$upper_shape <- rep(NaN, nrow(values))
  
  values$lower_mu <- rep(NaN, nrow(values))
  values$median_mu <- rep(NaN, nrow(values))
  values$upper_mu <- rep(NaN, nrow(values))
  for(i in 1:length(unique_temp)){
    values[i, "lower_rate"] <- quantile(fit_extract$rate_total_sporozoite[,i],0.025)
    values[i, "median_rate"] <- quantile(fit_extract$rate_total_sporozoite[,i],0.5)
    values[i, "upper_rate"] <- quantile(fit_extract$rate_total_sporozoite[,i],0.975)
    values[i, "lower_shape"] <- quantile(fit_extract$shape_total_sporozoite[,i],0.025)
    values[i, "median_shape"] <- quantile(fit_extract$shape_total_sporozoite[,i],0.5)
    values[i, "upper_shape"] <- quantile(fit_extract$shape_total_sporozoite[,i],0.975)
    values[i, "lower_mu"] <- quantile(fit_extract$mu_total_sporozoite[,i],0.025)
    values[i, "median_mu"] <- quantile(fit_extract$mu_total_sporozoite[,i],0.5)
    values[i, "upper_mu"] <- quantile(fit_extract$mu_total_sporozoite[,i],0.975)
  }
  return(values)
}



extract_rate_single <- function(temp_indexes, unique_temp_original){
  
  rate_single_values <- as.data.frame(unique_temp_original[temp_indexes])
  colnames(rate_single_values) <- c("temp")
  
  rate_single_values$lower_rate <- rep(NaN, nrow(rate_single_values))
  rate_single_values$median_rate <- rep(NaN, nrow(rate_single_values))
  rate_single_values$upper_rate <- rep(NaN, nrow(rate_single_values))
  
  rate_single_values$lower_shape <- rep(NaN, nrow(rate_single_values))
  rate_single_values$median_shape <- rep(NaN, nrow(rate_single_values))
  rate_single_values$upper_shape <- rep(NaN, nrow(rate_single_values))
  
  rate_single_values$lower_shape <- rep(NaN, nrow(rate_single_values))
  rate_single_values$median_shape <- rep(NaN, nrow(rate_single_values))
  rate_single_values$upper_shape <- rep(NaN, nrow(rate_single_values))
  
  for(i in 1:length(temp_indexes)){
    
    index <- temp_indexes[i]
    fit_single_rate <- readRDS(file = paste0("results/mSOS_single_temperature_model_fit_new_temp_",unique_temp_original[index]))
    fit_extract_single_rate <- rstan::extract(fit_single_rate)
    
    if(index == 3){ # for 27 degrees
      
      parameter_rate_single <- data.frame(rate_sporozoite = fit_extract_single_rate$rate_total_sporozoite,
                                          shape_sporozoite = fit_extract_single_rate$shape_total_sporozoite,
                                          mu_total = fit_extract_single_rate$mu_total_sporozoite)
      
    } else{
      parameter_rate_single <- data.frame(rate_sporozoite = fit_extract_single_rate$rate_sporozoite,
                                          shape_sporozoite = fit_extract_single_rate$shape_sporozoite)
      
      parameter_rate_single$mu_total <- parameter_rate_single[,"shape_sporozoite"] / parameter_rate_single[,"rate_sporozoite"]
    }
      
    rate_single_values[i,"lower_rate"] <- quantile(parameter_rate_single[,"rate_sporozoite"], 0.025)
    rate_single_values[i,"median_rate"] <- quantile(parameter_rate_single[,"rate_sporozoite"], 0.5)
    rate_single_values[i,"upper_rate"] <- quantile(parameter_rate_single[,"rate_sporozoite"], 0.975)
      
    rate_single_values[i,"lower_shape"] <- quantile(parameter_rate_single[,"shape_sporozoite"], 0.025)
    rate_single_values[i,"median_shape"] <- quantile(parameter_rate_single[,"shape_sporozoite"], 0.5)
    rate_single_values[i,"upper_shape"] <- quantile(parameter_rate_single[,"shape_sporozoite"], 0.975)
      
    rate_single_values[i,"lower_mu_total"] <- quantile(parameter_rate_single[,"mu_total"], 0.025)
    rate_single_values[i,"median_mu_total"] <- quantile(parameter_rate_single[,"mu_total"], 0.5)
    rate_single_values[i,"upper_mu_total"] <- quantile(parameter_rate_single[,"mu_total"], 0.975)
  
  }
  
  rate_single_values$label_a <- rep("A", nrow(rate_single_values))
  rate_single_values$label_b <- rep("B", nrow(rate_single_values))
  return(rate_single_values)
}

single_temp_model_param_plot <- function(to_go_rate_single_values, blank_theme){
  rate_G_to_S_plot <- ggplot() +
    geom_pointrange(data = to_go_rate_single_values, aes(x = temp, y = median_rate, ymin = lower_rate, 
                                                         ymax = upper_rate)) + blank_theme +
    ylab(parse(text = "beta[GS]")) + scale_x_continuous(limits = c(21,35), breaks = seq(21,35,4)) + 
    facet_wrap(.~label_b) +  theme(strip.text = element_text(face="bold", size=20), axis.title.y = element_text(size = 24))
  
  shape_G_to_S_plot <- ggplot() +
    geom_pointrange(data = to_go_rate_single_values, aes(x = temp, y = median_shape, ymin = lower_shape, 
                                                         ymax = upper_shape)) +  blank_theme +
    ylab(parse(text = "alpha[GS]")) + scale_x_continuous(limits = c(21,35), breaks = seq(21,35,4)) + 
    facet_wrap(.~label_a) +  theme(strip.text = element_text(face="bold", size=20), axis.title.y = element_text(size = 24))
  
  p <- plot_grid(shape_G_to_S_plot + theme(axis.title.x = element_blank()), 
                 rate_G_to_S_plot + theme(axis.title.x = element_blank())) 
  
  x_axis <- ggdraw() + draw_label("Temperature / °C", color = "black", fontface = "plain", size = 20,
                                  vjust = 0.25) + theme(plot.margin = margin(0, 0, 0, 0))
  
  plot_grid(p, x_axis, ncol = 1, rel_heights = c(1, 0.05))
}

#############################################################################
################ gamma distribution parameters at 27 ########################
#############################################################################

generate_PDF <- function(time, shape, rate){
  PDF <- data.frame("times" = seq(0,time, 0.1), "lower" = rep(NaN, length(seq(0,time, 0.1))), 
                      "median" = rep(NaN, length(seq(0,time, 0.1))), "upper" = rep(NaN, length(seq(0,time, 0.1))))
  
  for(i in 1:nrow(PDF)){
    placeholder <- dgamma(PDF[i,"times"], shape, rate)
    PDF[i, "lower"] = quantile(placeholder, 0.025)[[1]]
    PDF[i, "median"] = quantile(placeholder, 0.5)[[1]]
    PDF[i, "upper"] = quantile(placeholder, 0.975)[[1]]
    rm(list = c("placeholder"))
  }
  return(PDF)
}

generate_CDF <- function(time, shape, rate){
  PDF <- data.frame("times" = seq(0,time, 0.1), "lower" = rep(NaN, length(seq(0,time, 0.1))), 
                    "median" = rep(NaN, length(seq(0,time, 0.1))), "upper" = rep(NaN, length(seq(0,time, 0.1))))
  
  for(i in 1:nrow(PDF)){
    placeholder <- pgamma(PDF[i,"times"], shape, rate)
    PDF[i, "lower"] = quantile(placeholder, 0.025)[[1]]
    PDF[i, "median"] = quantile(placeholder, 0.5)[[1]]
    PDF[i, "upper"] = quantile(placeholder, 0.975)[[1]]
    rm(list = c("placeholder"))
  }
  return(PDF)
}

generate_density_df <- function(fit_27_extract){
  PDF_O <- generate_PDF(10, fit_27_extract$shape_oocyst, fit_27_extract$rate_oocyst)
  PDF_S <- generate_PDF(20, fit_27_extract$shape_sporozoite, fit_27_extract$rate_sporozoite)
  PDF_T <- generate_PDF(25, fit_27_extract$shape_total_sporozoite, fit_27_extract$rate_total_sporozoite)
  PDF_O$stage <- rep("G to O", nrow(PDF_O))
  PDF_S$stage <- rep("O to S", nrow(PDF_S))
  PDF_T$stage <- rep("G to S", nrow(PDF_T))
  PDF <- rbind(PDF_O, PDF_S, PDF_T)
  PDF$stage <- factor(PDF$stage, levels = c("G to O", "O to S", "G to S"))
  PDF$density <- rep("A (PDF)", nrow(PDF))

  CDF_O <- generate_CDF(10, fit_27_extract$shape_oocyst, fit_27_extract$rate_oocyst)
  CDF_S <- generate_CDF(20, fit_27_extract$shape_sporozoite, fit_27_extract$rate_sporozoite)
  CDF_T <- generate_CDF(25, fit_27_extract$shape_total_sporozoite, fit_27_extract$rate_total_sporozoite)
  CDF_O$stage <- rep("G to O", nrow(CDF_O))
  CDF_S$stage <- rep("O to S", nrow(CDF_S))
  CDF_T$stage <- rep("G to S", nrow(CDF_T))
  CDF <- rbind(CDF_O, CDF_S, CDF_T)
  CDF$stage <- factor(CDF$stage, levels = c("G to O", "O to S", "G to S"))
  CDF$density <- rep("B (CDF)", nrow(CDF))
  
  prob_df <- rbind(PDF, CDF)
}

prob_plot <- function(prob_df){
  ggplot(data = prob_df) + 
  geom_ribbon(aes(x = times, ymin = lower, ymax = upper, group = stage), fill = "grey55", alpha = 0.4) +
  geom_line(aes(x = times, y = median, group = stage, colour = stage), size = 0.75) +
  xlab("Days") + ylab("Density") + 
  scale_colour_manual(values = c("#E69F00", "#56B4E9", "#CC79A7")) +
  blank_theme + facet_wrap(.~density, scales = "free") + theme(strip.text = element_text(face="bold", size=20)) +
  labs(group = "Life stage \n transition", colour = "Life stage \n transition")
}


