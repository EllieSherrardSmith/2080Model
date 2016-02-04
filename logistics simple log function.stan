data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real<lower=0, upper=1> alpha;
  real<lower=0> eps;
}
model {
  for (n in 1:N)
    y[n] ~ normal((1 / (1 + exp(alpha * x[n]))), eps); 
    increment_log_prob(-log(eps));     //log prior for p(eps) proportion to 1/eps
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    log_lik[n] <- normal_log(y[n], (1 / (1 + exp(alpha * x[n]))),eps);
  }
}