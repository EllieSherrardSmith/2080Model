data {
  int<lower=0> N;
  vector[N] x1;
  vector[N] x2;
  vector[N] x3;
  vector[N] y;
}
parameters {
  real<lower=0.7, upper=0.9> alpha;
  real beta1;
  real beta2;
  real beta3;

  real<lower=0> sigma;
}
model {
  y ~ normal(alpha + beta1 * x1 + beta2 * x2 + beta3 * x3, sigma);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    log_lik[n] <- normal_log(y[n], (alpha + beta1 * x1[n] + beta2 * x2[n] + beta3 * x3[n]), sigma);
  }
}