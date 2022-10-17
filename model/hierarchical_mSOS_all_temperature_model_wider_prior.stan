// mSOS: a multiscale stochastic model of sporogony - all temperature model
// Author: Isaac J Stopard
// Version: 0.01 
// Last updated: 29/04/2020
// Notes: within vector model of malaria parasite population dynamics
//        fits the malaria parasite dissection data for all temperatures simultaneously
//        fits to oocyst prevalence, oocyst intensity (the parasite load within individual mosquitoes), sporozoite prevalence and survival data simultaneously
//        zero truncated negative binomial distribution is used to model the distribution of parasites within the mosquitoes
//        parasite development times are modelled as gamma distributed

functions{
  
  // incomplete regularised gamma function
  real regularized_gamma(real a, real z_0, real z_1){
    real result = gamma_p(a, z_1) - gamma_p(a, z_0);
    return result;
  }
  
  // Parasite prevalence estimate (oocyst or sporozoite) 
  // inputs: t (time), mu_NB (negative binomial distribution mean parameter), k (negative binomial distribution overdispersion parameter), 
  //             a (gamma distribution shape parameter) and b (gamma distribution rate parameter)
  real parasite_prevalence(real t, real mu_NB, real k, real a, real b){
    real result;
    result = (1 - (k^k) * ((k + mu_NB * regularized_gamma(a, 0, t * b))^-k)) / (1 - (k / (k + mu_NB))^k);
    return result;
  
  }
  
  // log PDF for the oocyst intensity data
  // inputs: intensity (observed oocyst load within individual mosquitoes), t (time), mu_NB (negative binomial distribution mean parameter),
  //         k (negative binomial distribution overdispersion parameter), a (gamma distribution shape parameter), b (gamma distribution rate parameter),
  //         delta (probability of viable infection)
  real parasite_intensity_PDF(real intensity, real t, real mu_NB, real k, real a, real b, real delta){
    real result;
    real theta = regularized_gamma(a, 0, t * b);
    if(intensity == 0){
      result = log((1 - delta) + delta * ((k / (k + mu_NB))^k - (k / (k + theta * mu_NB))^k) / (-1 + (k / (k + mu_NB))^k));
    } else{
      result = log(delta * ((k^k * (theta * mu_NB)^intensity * (k + theta * mu_NB)^(-intensity - k) * exp(lchoose(-1 + intensity + k, -1 + k))) / (1 - (k / (k + mu_NB))^k)));
      }
    return result;
  }
  
  // functions to describe the relationship between temperature and different model parameters
  real logistic_function(real x){
    real result = 1 / (1 + exp(-x));
    return result;
  }
  
  real quadratic_function(real a, real b, real c, real temp){
    real result = (a * temp^2) + (b * temp) + c;
    return result;
  }
  
  real linear_function(real m, real c, real temp, real error){
    real result = (m * temp) + c + error;
    return result;
  }

  // survival functions
  real cox(real intercept, real beta_inf, real beta_temp, real infec, real temp, real error){
    // real result = exp(intercept + (beta_inf * infec) + (beta_temp * temp)); not log
    real result = intercept + (beta_inf * infec) + (beta_temp * temp) + error;
    return result;
  }
  
  real baseline_hazard(real t, real a, real b){
    // real result = a * exp(b * t); // not log
    real result = log(a) + (b * t);
    return result;
  }
  
  real S_t(real t, real a, real b, real cox_){
    real result = (a/b)*(1-exp(b*t))*exp(cox_);
    return result;
  }
  
  real log_likelihood_survival(real delta, int infection_status, int censored, real base_hazard, real cox_u, real cox_i, real S_t_i, real S_t_u){
    real log_likelihood;
    //if(censored == 1){
    //    if(infection_status == 1){
    //      log_likelihood = log((delta * base_hazard * cox_i * S_t_i) + 
    //                                      ((1-delta) * base_hazard * cox_u * S_t_u));
    //    } else if(infection_status < 1){
    //      log_likelihood = log(base_hazard * cox_u * S_t_u);
    //    }
    //} else if(censored < 1){
    //    if(infection_status == 1){
    //        log_likelihood = log((delta * S_t_i) + ((1-delta) * S_t_u)); 
    //  } else if(infection_status < 1){
    //        log_likelihood = log(S_t_u);
    //  }
    //}
    
      if(censored == 1){
        if(infection_status == 1){
          log_likelihood = log_sum_exp((log(delta) + base_hazard + cox_i + S_t_i),(log(1-delta) + base_hazard + cox_u + S_t_u));
        } else if(infection_status < 1){
          log_likelihood = base_hazard + cox_u + S_t_u;
        }
    } else if(censored < 1){
        if(infection_status == 1){
          log_likelihood = log_sum_exp((log(delta) + S_t_i), (log(1-delta) + S_t_u));
      } else if(infection_status < 1){
           log_likelihood = S_t_u;
      }
    }
    return log_likelihood;
  }
}

data{
  int length_oocysts;
  int oocyst_total_sampled[length_oocysts]; // n
  int oocyst_total_positive[length_oocysts]; // X
  int oocyst_index_temp[length_oocysts]; // temperature index
  int oocyst_exp[length_oocysts];
  real oocyst_time[length_oocysts]; // temperature index and day post infection index
  
  int length_sporozoites;
  int sporozoite_total_sampled[length_sporozoites]; // n
  int sporozoite_total_positive[length_sporozoites]; // X
  int sporozoite_index_temp[length_sporozoites];
  int sporozoite_exp[length_sporozoites];
  real sporozoite_time[length_sporozoites];
  
  int length_unique_temp;
  real unique_temp[length_unique_temp];
  
  int length_unique_exp;
  
  int length_ppd_times;
  vector[length_ppd_times] PPD_times;
  
  // survival data
  int length_survival_data;
  int survival_index[length_survival_data];
  int length_unique_survival_indices;
  vector[length_unique_survival_indices] unique_survival_times;
  int unique_survival_temp_index[length_unique_survival_indices];
  int unique_survival_exp[length_unique_survival_indices];
  int unique_survival_infection_status[length_unique_survival_indices];
  int unique_survival_censored[length_unique_survival_indices];
  
  // oocyst intensity data
  int length_oocyst_intensity;
  int oocyst_intensity_index[length_oocyst_intensity];
  
  int length_unique_oocyst_intensity_indices;
  real unique_oocyst_intensity_time[length_unique_oocyst_intensity_indices];
  real unique_oocyst_intensity[length_unique_oocyst_intensity_indices];
  int unique_oocyst_intensity_temp_index[length_unique_oocyst_intensity_indices];
  int unique_oocyst_intensity_exp[length_unique_oocyst_intensity_indices];
  real intercept;
}

parameters{
  // development time parameters
  real<lower=0> shape_oocyst;
  real m_rate_oocyst;
  real<lower=0> c_rate_oocyst;
  
  real<lower=0>  shape_sporozoite;
  real<lower=0> rate_sporozoite;
  
  // parasite intensity population parameters
  real<lower=0> mu_NB;
  real<lower=0> k_NB;
  
  // proportion of mosquitoes infected parameters (probability of viable infection)
  real m_delta;
  real<lower=0> c_delta;
  
  // mosquito survival parameters
  real<lower=0, upper=5> a;
  real<lower=0, upper=5> b;
  real<lower=-2.5, upper=2.5> beta_inf;
  real<lower=-2.5, upper=2.5> beta_temp;
  // real<lower=-10, upper=5> intercept;
  
      // errors
  real error_delta[length_unique_exp];
  real error_survival[length_unique_exp];
  real<lower=0> sigma_error_delta;
  real<lower=0> sigma_error_survival;
}

transformed parameters{
  
  real<lower=0, upper=1> delta[length_unique_temp, length_unique_exp];
  real<lower=0> rate_oocyst[length_unique_temp];
  
  real<lower=0> mu_total_sporozoite[length_unique_temp];
  real<lower=0> sigma_sq_sporozoite[length_unique_temp];
  
  real<lower=0> shape_total_sporozoite[length_unique_temp];
  real<lower=0> rate_total_sporozoite[length_unique_temp];
  
  vector<lower=0>[length_oocysts] RR_oocyst;
  vector<lower=0>[length_sporozoites] RR_sporozoite;
  
  // oocyst intensity data
  vector<lower=0>[length_unique_oocyst_intensity_indices] RR_oocyst_intensity_unique;
  vector[length_unique_oocyst_intensity_indices] oocyst_intensity_probability_lookup;
  
  // survival data
  vector[length_unique_survival_indices] baseline_hazard_lookup;
  vector[length_unique_survival_indices] cox_infected_lookup;
  vector[length_unique_survival_indices] cox_uninfected_lookup;
  vector[length_unique_survival_indices] S_t_infected_lookup;
  vector[length_unique_survival_indices] S_t_uninfected_lookup;
  vector[length_unique_survival_indices] survival_likelihood_lookup;
  
  
  vector[length_oocysts] oocyst_likelihoods;
  vector[length_oocysts] oocyst_prevalence;
  
  vector[length_sporozoites] sporozoite_likelihoods;
  vector[length_sporozoites] sporozoite_prevalence;
  
  vector[length_oocyst_intensity] oocyst_intensity_likelihoods;
  
  // survival analysis
  real log_likelihood_survival_data[length_survival_data];
  
  for(i in 1:length_unique_temp){
    for(j in 1:length_unique_exp){
    delta[i,j] = logistic_function(linear_function(m_delta, c_delta, unique_temp[i], error_delta[j]));
    }
  }
  
  for(i in 1:length_unique_temp){
    rate_oocyst[i] = linear_function(m_rate_oocyst, c_rate_oocyst, unique_temp[i], 0); // no error is included for differences in parasite development rate
    mu_total_sporozoite[i] = (rate_oocyst[i] * shape_sporozoite + rate_sporozoite * shape_oocyst) / (rate_oocyst[i] * rate_sporozoite);
    sigma_sq_sporozoite[i] = (rate_oocyst[i]^2 * shape_sporozoite + rate_sporozoite^2 * shape_oocyst) / (rate_oocyst[i]^2 * rate_sporozoite^2);
    shape_total_sporozoite[i] = mu_total_sporozoite[i]^2 / sigma_sq_sporozoite[i];
    rate_total_sporozoite[i] = mu_total_sporozoite[i] / sigma_sq_sporozoite[i];
  }
  
  for(i in 1:length_oocysts){
  RR_oocyst[i] = exp(S_t(oocyst_time[i], a, b, cox(intercept, beta_inf, beta_temp, 1.0, unique_temp[oocyst_index_temp[i]], error_survival[oocyst_exp[i]]))) / 
                    exp(S_t(oocyst_time[i], a, b, cox(intercept, beta_inf, beta_temp, 0, unique_temp[oocyst_index_temp[i]], error_survival[oocyst_exp[i]])));
  }
  
  for(i in 1:length_sporozoites){
    RR_sporozoite[i] = exp(S_t(sporozoite_time[i], a, b, cox(intercept, beta_inf, beta_temp, 1.0, unique_temp[sporozoite_index_temp[i]], error_survival[sporozoite_exp[i]]))) / 
                       exp(S_t(sporozoite_time[i], a, b, cox(intercept, beta_inf, beta_temp, 0, unique_temp[sporozoite_index_temp[i]], error_survival[sporozoite_exp[i]])));
  }
  
  for(i in 1:length_unique_oocyst_intensity_indices){
    RR_oocyst_intensity_unique[i] = exp(S_t(unique_oocyst_intensity_time[i], a, b, cox(intercept, beta_inf, beta_temp, 1.0, unique_temp[unique_oocyst_intensity_temp_index[i]], error_survival[unique_oocyst_intensity_exp[i]]))) / 
                    exp(S_t(unique_oocyst_intensity_time[i], a, b, cox(intercept, beta_inf, beta_temp, 0, unique_temp[unique_oocyst_intensity_temp_index[i]], error_survival[unique_oocyst_intensity_exp[i]])));
  }
  
  // survival data
  for(i in 1:length_unique_survival_indices){
    baseline_hazard_lookup[i] = baseline_hazard(unique_survival_times[i], a, b);
    cox_infected_lookup[i] = cox(intercept, beta_inf, beta_temp, 1, unique_temp[unique_survival_temp_index[i]], error_survival[unique_survival_exp[i]]);
    cox_uninfected_lookup[i] = cox(intercept, beta_inf, beta_temp, 0, unique_temp[unique_survival_temp_index[i]], error_survival[unique_survival_exp[i]]);
    S_t_infected_lookup[i] = S_t(unique_survival_times[i], a, b, cox_infected_lookup[i]);
    S_t_uninfected_lookup[i] = S_t(unique_survival_times[i], a, b, cox_uninfected_lookup[i]);
    survival_likelihood_lookup[i] = log_likelihood_survival(delta[unique_survival_temp_index[i], unique_survival_exp[i]], unique_survival_infection_status[i], 
                                    unique_survival_censored[i], baseline_hazard_lookup[i], cox_uninfected_lookup[i], 
                                    cox_infected_lookup[i], S_t_infected_lookup[i], S_t_uninfected_lookup[i]);
  }
  
  // oocyst intensity data
  for(i in 1:length_unique_oocyst_intensity_indices){
    oocyst_intensity_probability_lookup[i] = parasite_intensity_PDF(unique_oocyst_intensity[i], unique_oocyst_intensity_time[i], 
                                            mu_NB, k_NB, shape_oocyst, rate_oocyst[unique_oocyst_intensity_temp_index[i]], 
                                            (delta[unique_oocyst_intensity_temp_index[i], unique_oocyst_intensity_exp[i]] * RR_oocyst_intensity_unique[i]));
  }
  
  // calculations
  for(i in 1:length_oocysts){
    oocyst_prevalence[i] = parasite_prevalence(oocyst_time[i], mu_NB, k_NB, shape_oocyst, rate_oocyst[oocyst_index_temp[i]]);
    oocyst_likelihoods[i] = binomial_lpmf(oocyst_total_positive[i] | oocyst_total_sampled[i], (delta[oocyst_index_temp[i], oocyst_exp[i]] * RR_oocyst[i] * oocyst_prevalence[i]));
  }
  
  for(i in 1:length_sporozoites){
    sporozoite_prevalence[i] = parasite_prevalence(sporozoite_time[i], mu_NB, k_NB, shape_total_sporozoite[sporozoite_index_temp[i]], rate_total_sporozoite[sporozoite_index_temp[i]]);
    sporozoite_likelihoods[i] = binomial_lpmf(sporozoite_total_positive[i] | sporozoite_total_sampled[i], (delta[sporozoite_index_temp[i], sporozoite_exp[i]] * RR_sporozoite[i] * sporozoite_prevalence[i]));
  }
  
  // survival analysis
  for(i in 1:length_survival_data){
    log_likelihood_survival_data[i] = survival_likelihood_lookup[survival_index[i]];
  }
  
  for(i in 1:length_oocyst_intensity){
    oocyst_intensity_likelihoods[i] = oocyst_intensity_probability_lookup[oocyst_intensity_index[i]];
  }
  
}

model{
  // priors
  shape_oocyst ~ normal(15, 4.5);
  m_rate_oocyst ~ normal(0, 2.5);
  c_rate_oocyst ~ normal(3.0, 2.5);
  
  shape_sporozoite ~ normal(15, 4.5); 
  rate_sporozoite ~ normal(1.875, 2.5);
  
  m_delta ~ normal(0, 2.5);
  c_delta ~ normal(0.8, 2.5); 
  
  mu_NB ~ normal(3.0, 2.5);
  k_NB ~  normal(0.1, 2.5) ;
  
  // survival analysis
  a ~ normal(0.05, 2.5);
  b ~ normal(0.05, 2.5);
  beta_inf ~ normal(0.5, 1.0);
  beta_temp ~ normal(0.5, 1.0);
  //intercept ~ normal(-3, 2.5);
  
  for(i in 1:length_unique_exp){
    error_delta[i] ~ normal(0, sigma_error_delta);
    error_survival[i] ~ normal(0, sigma_error_survival);
  }
  
  sigma_error_delta ~ normal(0, 0.25);
  sigma_error_delta ~ normal(0, 0.25);
  
  target += (sum(oocyst_likelihoods) + sum(sporozoite_likelihoods) + sum(log_likelihood_survival_data) + sum(oocyst_intensity_likelihoods));
}

generated quantities{
  matrix[length_ppd_times, length_unique_temp] RR_ppd;
  matrix[length_ppd_times, length_unique_temp] oocyst_prevalence_ppd;
  matrix[length_ppd_times, length_unique_temp] sporozoite_prevalence_ppd;
  matrix[length_ppd_times, length_unique_temp] St_uninfected_ppd;
  matrix[length_ppd_times, length_unique_temp] St_infected_ppd;
  matrix[length_ppd_times, length_unique_temp] St_infected_blood_fed_ppd;
  matrix[length_ppd_times, length_unique_temp] oocyst_intensity_ppd;
  real mean_oocyst = mu_NB / (1 - (k_NB / (k_NB + mu_NB))^k_NB);
  
  real<lower=0,upper=1> delta_no_error[length_unique_temp];
  
  for(i in 1:length_unique_temp){
    delta_no_error[i] = logistic_function(linear_function(m_delta, c_delta, unique_temp[i], 0)); // no error parameter
  }
  
  for(i in 1:length_ppd_times){
    for(j in 1:length_unique_temp){
      St_uninfected_ppd[i,j] = exp(S_t(PPD_times[i], a, b, cox(intercept, beta_inf, beta_temp, 0, unique_temp[j], 0))); // no error parameter
      St_infected_ppd[i,j] = exp(S_t(PPD_times[i], a, b, cox(intercept, beta_inf, beta_temp, 1.0, unique_temp[j], 0))); // no error parameter
      St_infected_blood_fed_ppd[i,j] = (exp(S_t(PPD_times[i], a, b, cox(intercept, beta_inf, beta_temp, 1.0, unique_temp[j], 0))) * delta_no_error[j]) +
                                        (exp(S_t(PPD_times[i], a, b, cox(intercept, beta_inf, beta_temp, 0, unique_temp[j], 0))) * (1-delta_no_error[j])); // no error parameter
      RR_ppd[i,j] = St_infected_ppd[i,j]/St_uninfected_ppd[i,j]; // no error parameter
      oocyst_prevalence_ppd[i,j] = parasite_prevalence(PPD_times[i], mu_NB, k_NB, shape_oocyst, rate_oocyst[j]) * delta_no_error[j] * RR_ppd[i,j];
      sporozoite_prevalence_ppd[i,j] = parasite_prevalence(PPD_times[i], mu_NB, k_NB, shape_total_sporozoite[j], rate_total_sporozoite[j]) * delta_no_error[j] * RR_ppd[i,j];
      oocyst_intensity_ppd[i,j] = delta_no_error[j] * RR_ppd[i, j] * mean_oocyst * regularized_gamma(shape_oocyst, 0, (PPD_times[i] * rate_oocyst[j]));
    }
  }
}

