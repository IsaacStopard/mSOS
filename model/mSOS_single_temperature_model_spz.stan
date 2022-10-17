// mSOS: a multiscale stochastic model of sporogony - single temperature model
// Author: Isaac J Stopard
// Version: 0.01 
// Last updated: 29/04/2020
// Notes: within vector model of malaria parasite population dynamics
//        fits the malaria parasite dissection data for single temperatures
//        fits to sporozoite prevalence and survival data simultaneously
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
  
  // log PDF for the parasite prevalence data
  // inputs: n (total mosquito sample), X (number of parasite positive mosquitoes), delta (probability of viable infection), 
  //         prevalence (model parasite prevalence output)
  real compound_binomial_PDF(int n, int X, real delta, real prevalence){
    real result;
    
    // calculated not on the log scale
    //if(X == 0){ // none showing parasites
    //  result = (1 - (prevalence * delta)) ^ n;
    //} else if(n == X){ // all showing parasites
    //  result = ((1 - prevalence) ^ (n - X)) * (prevalence ^ X) * (delta ^ n) * exp(lchoose(n, X));
    //} else if(n > X && X > 0){
    //  result = ((1 - delta) ^ (n - X)) * ((prevalence * delta) ^ X) * (((-1 + (prevalence * delta))/(-1 + delta)) ^ (n - X)) * exp(lchoose(n, X));
    //}
    //result = log(result); // PDF must be on the log scale
    
    // calculated on the log scale
    if(X == 0){ // none showing parasites
      result = n * log(1 - (prevalence * delta));
    } else if(n == X){ // all showing parasites
      result = ((n - X) * log(1 - prevalence)) + (X * log(prevalence)) + (n * log(delta)) + lchoose(n, X);
    } else if(n > X && X > 0){
      result = ((n - X) * log(1 - delta)) + (X * log(prevalence * delta)) + ((n - X) * log((-1 + (prevalence * delta))/(-1 + delta))) + lchoose(n, X);
    }
    return result;
  }

  // survival functions
  real cox(real intercept, real beta_inf, real infec){
    real result = intercept + (beta_inf * infec);
    return result;
  }
  
  real baseline_hazard(real t, real a, real b){
    // real result = a * exp(b * t);
    real result = log(a) + (b * t);
    return result;
  }
  
  real S_t(real t, real a, real b, real cox_){
    //real result = exp((a/b)*(1-exp(b*t)))^cox_;
    real result = (a/b)*(1-exp(b*t))*exp(cox_);
    return result;
  }
  
  real log_likelihood_survival(real delta, int infection_status, int censored, real base_hazard, real cox_u, real cox_i, real S_t_i, real S_t_u){
    real log_likelihood;
    //if(censored == 1){
    //    if(infection_status == 1){
    //      log_likelihood = log((delta * base_hazard * S_t_i * cox_i) + 
    //                                      ((1-delta) * base_hazard * S_t_u * cox_u));
    //    } else if(infection_status < 1){
    //      log_likelihood = log(base_hazard * S_t_u * cox_u);
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
  int length_sporozoites;
  int sporozoite_total_sampled[length_sporozoites]; // n
  int sporozoite_total_positive[length_sporozoites]; // X
  real sporozoite_time[length_sporozoites];
  
  int length_ppd_times;
  vector[length_ppd_times] PPD_times;
  
  // survival data
  int length_survival_data;
  int survival_index[length_survival_data];
  int length_unique_survival_indices;
  vector[length_unique_survival_indices] unique_survival_times;
  int unique_survival_infection_status[length_unique_survival_indices];
  int unique_survival_censored[length_unique_survival_indices];
  real intercept;
}

parameters{
  // development time parameters
  real<lower=0>  shape_sporozoite;
  real<lower=0> rate_sporozoite;
  
  // parasite intensity population parameters
  real<lower=0 > mu_NB;
  real<lower=0 > k_NB; // real<lower=0, upper=50> k_NB; // for stephensi only data
  
  // proportion of mosquitoes infected parameters
  // real<lower=0, upper=1> delta_base;
  // real a_delta;
  real<lower=0, upper=1> delta;

  // mosquito survival parameters
  real<lower=0, upper=5> a;
  real<lower=0, upper=5> b;
  real<lower=-5, upper=5> beta_inf;
  //real<lower=-10, upper=5> intercept;
}

transformed parameters{
  
  vector<lower=0>[length_sporozoites] RR_sporozoite;

  // survival data
  vector[length_unique_survival_indices] baseline_hazard_lookup;
  vector[length_unique_survival_indices] cox_infected_lookup;
  vector[length_unique_survival_indices] cox_uninfected_lookup;
  vector[length_unique_survival_indices] S_t_infected_lookup;
  vector[length_unique_survival_indices] S_t_uninfected_lookup;
  vector[length_unique_survival_indices] survival_likelihood_lookup;

  for(i in 1:length_sporozoites){
    RR_sporozoite[i] = exp(S_t(sporozoite_time[i], a, b, cox(intercept, beta_inf, 1.0))) / 
                       exp(S_t(sporozoite_time[i], a, b, cox(intercept, beta_inf, 0)));
  }
  
  // survival data
  for(i in 1:length_unique_survival_indices){
    baseline_hazard_lookup[i] = baseline_hazard(unique_survival_times[i], a, b);
    cox_infected_lookup[i] = cox(intercept, beta_inf, 1);
    cox_uninfected_lookup[i] = cox(intercept, beta_inf, 0);
    S_t_infected_lookup[i] = S_t(unique_survival_times[i], a, b, cox_infected_lookup[i]);
    S_t_uninfected_lookup[i] = S_t(unique_survival_times[i], a, b, cox_uninfected_lookup[i]);
    survival_likelihood_lookup[i] = log_likelihood_survival(delta, unique_survival_infection_status[i], 
                                    unique_survival_censored[i], baseline_hazard_lookup[i], cox_uninfected_lookup[i], 
                                    cox_infected_lookup[i], S_t_infected_lookup[i], S_t_uninfected_lookup[i]);
  }
}

model{

  vector[length_sporozoites] sporozoite_likelihoods;
  vector[length_sporozoites] sporozoite_prevalence;
  
  // survival analysis
  real log_likelihood_survival_data[length_survival_data];
  
  for(i in 1:length_sporozoites){
    sporozoite_prevalence[i] = parasite_prevalence(sporozoite_time[i], mu_NB, k_NB, shape_sporozoite, rate_sporozoite);
    // sporozoite_likelihoods[i] = compound_binomial_PDF(sporozoite_total_sampled[i], sporozoite_total_positive[i], delta * RR_sporozoite[i], sporozoite_prevalence[i]);
    sporozoite_likelihoods[i] = binomial_lpmf(sporozoite_total_positive[i] | sporozoite_total_sampled[i], (delta * RR_sporozoite[i] * sporozoite_prevalence[i]));
  }
  
  // survival analysis
  for(i in 1:length_survival_data){
    log_likelihood_survival_data[i] = survival_likelihood_lookup[survival_index[i]];
  }
  
  // priors
  shape_sporozoite ~ normal(15, 2.5);
  rate_sporozoite ~ normal(1.875, 2.5); // normal(2, 2.5) // for stephensi only data
  
  //delta_base ~ beta(5.5, 2.5);
  delta ~ beta(5.5, 2.5);
 
  mu_NB ~ normal(3.0, 2.5); // normal(3, 2.5); // for stephensi only data
  k_NB ~  normal(0.1, 2.5) ; // normal(5, 2.5) for stephensi only data
  
  // survival analysis
  a ~ normal(0.05, 2.5);
  b ~ normal(0.05, 2.5);
  beta_inf ~ normal(0.5, 1.0);
  //intercept ~ normal(-3, 2.5);
  
  target += sum(sporozoite_likelihoods) + sum(log_likelihood_survival_data);
}

generated quantities{
  vector[length_ppd_times] RR_ppd;
  vector[length_ppd_times] sporozoite_prevalence_ppd;
  vector[length_ppd_times] St_uninfected_ppd;
  vector[length_ppd_times] St_infected_ppd;
  vector[length_ppd_times] St_infected_blood_fed_ppd;
  
  for(i in 1:length_ppd_times){
      St_uninfected_ppd[i] = exp(S_t(PPD_times[i], a, b, cox(intercept, beta_inf, 0)));
      St_infected_ppd[i] = exp(S_t(PPD_times[i], a, b, cox(intercept, beta_inf, 1.0)));
      St_infected_blood_fed_ppd[i] = (exp(S_t(PPD_times[i], a, b, cox(intercept, beta_inf, 1.0))) * delta) + 
                                      (exp(S_t(PPD_times[i], a, b, cox(intercept, beta_inf, 0))) * (1 - delta));
      RR_ppd[i] = St_infected_ppd[i]/St_uninfected_ppd[i];
      sporozoite_prevalence_ppd[i] = parasite_prevalence(PPD_times[i], mu_NB, k_NB, shape_sporozoite, rate_sporozoite) * delta *
                                       RR_ppd[i];
  }
}


