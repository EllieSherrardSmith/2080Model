data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real<lower=0> alpha;
  real<lower=0>  beta;
  real<lower=0> sigma;
}
model {
  for (n in 1:N)
    y[n] ~ normal(1 / (1 + exp(-(-alpha + beta x[n]))),sigma);
## not normal distribution (between 0 and 1)
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    log_lik[n] <- normal_log(y[n], 1 / (1 + exp(-(-alpha + beta x[n]))), sigma);
  }
}
1/(1+exp(-(-4.0777+1.5046* Hours)))