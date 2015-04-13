data{

  int<lower=0> N_st; ##number of treatments investigated (e.g. no treatment, remove 5%, 10% 20% etc)
  int<lower=0> N_counts; ##number of data counts "hosts with intensity of parasites" i.e. 21
  
  int<lower=0> para_count[N_counts,N_st]; ##the raw data counts of the oocysts for each group for controls
  int<lower=0> prev[N_counts,N_st];#the counts of individuals who are infected and those that are not
  
}
transformed data{
  int N_treatments;
  real offset;
  N_treatments <- 21;
  offset <- 2;
}
parameters{
  real rho__log_mean_para_count;
  real<lower=0> tau_Tx__log_mean_para_count;
  real rho_Tx__log_mean_para_count[N_treatments];
  
  real rho__log_od_para_count;
  real<lower=0> tau_Tx__log_od_para_count;
  real rho_Tx__log_od_para_count[N_treatments];

  ## Slopes and intercept of GLM modeling infection probability from
  ## parasite count distribution.  Common parameters for 
  ## different treatment.
  vector[2] beta_theta;
  real alpha_theta;
}
model{
  ## Log mean and log precision of parasite count NegBin distribution
  vector[N_st] log_mu_para_count;
  vector[N_st] log_phi_para_count;
  
  ## Logit infection probability for each control/treatment group
  vector[N_st] logit_theta;

  ## Priors
  rho__log_mean_para_count ~ normal(0, 5);
  tau_Tx__log_mean_para_count ~ cauchy(0, 2.5);
  rho_Tx__log_mean_para_count ~ normal(0, 1);
  
  rho__log_od_para_count ~ normal(0, 5);
  tau_Tx__log_od_para_count ~ cauchy(0, 2.5);
  rho_Tx__log_od_para_count ~ normal(0, 1);

  beta_theta ~ normal(0, 4);

  alpha_theta ~ normal(0, 4);


  ## Parasite counts
  for (g in 1:N_st) {
    log_mu_para_count[g] <-   rho__log_mean_para_count
                       + tau_Tx__log_mean_para_count
                         * rho_Tx__log_mean_para_count[g];
    log_phi_para_count[g] <-  rho__log_od_para_count
                       + tau_Tx__log_od_para_count
                         * rho_Tx__log_od_para_count[g];
  }
  
  for(n in 1:N_counts)
      para_count[n] ~ neg_binomial_2_log(log_mu_para_count, exp(log_phi_para_count));
  

  ## GLM of infection probability in place of intractable convolution
  for (n in 1:N_st)
    logit_theta[n] <-    beta_theta[1] * log_mu_para_count[n]
                         + beta_theta[2] * log_phi_para_count[n]
                         + alpha_theta;
  
  ## Host infection prevalence measurements
  for(n in 1:N_st) {
    prev[n] ~ bernoulli_logit(logit_theta);

  }
}
generated quantities {
  int<lower=0> sim_para_count[N_st];
  int<lower=0> sim_prev[N_st];
  real<lower=0, upper=1> theta[N_st];
  
  {
  for (g in 1:N_st) {
      real log_mu_para_count;
      real log_phi_para_count;
      
      // Simulate oocyst count measurement
      log_mu_para_count <-  rho__log_mean_para_count
                   + tau_Tx__log_mean_para_count
                   * rho_Tx__log_mean_para_count[g];
      log_phi_para_count <-  rho__log_od_para_count
                    + tau_Tx__log_od_para_count
                    * rho_Tx__log_od_para_count[g];
      
      sim_para_count[g]
        <- neg_binomial_2_log_rng(log_mu_para_count, exp(log_phi_para_count));
      
      theta[g] <- inv_logit(  beta_theta[1] * log_mu_para_count
                              + beta_theta[2] * log_phi_para_count
                              + alpha_theta);
  
    ## Simulate infection prevelance measurement
      sim_prev[g] <- binomial_rng(N_st, theta[g]);
      
    }
}
}