// logistic function
data {
  int<lower=1> N;    // rows of data
  vector[N] x;       // predictor
  real<lower=0> y[N]; // response
}

parameters {
  real phi; // 
    real beta;  // 
    real alpha;
}
model {
  // priors:
    phi ~ cauchy(0, 2);
  beta ~ normal(0, 2);
  alpha ~ normal(0, 2);
  
  // data model:
    for (n in 1:N)
      // y[n] ~ normal( (1 - ((alpha * x[n]) / sqrt(1 + pow(x[n], 1/beta)))), phi);
  y[n] ~ normal(((exp(alpha + beta * x[n])) / (1 + exp(alpha + beta * x[n])) ), phi);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    //    log_lik[n] <- normal_log(y[n], (1 - ((alpha * x[n]) / sqrt(1 + pow(x[n], 1/beta)))), phi);
    log_lik[n] <- normal_log(y[n], ((exp(alpha + beta * x[n])) / (1 + exp(alpha + beta * x[n])) ), phi);
  }
}