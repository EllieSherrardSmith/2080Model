data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real<lower=0.1, upper=0.5> alpha;
  real<lower=0.1, upper=0.9>  beta;
  real<lower=0, upper=1> delta;
  real<lower=0> sigma;
}
model {
  for (n in 1:N)
    y[n] ~ normal((alpha/beta) * exp(-exp(delta - beta * x[n])),sigma);
## not normal distribution (between 0 and 1)
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    log_lik[n] <- normal_log(y[n], (alpha/beta) * exp(-exp(delta - beta * x[n])), sigma);
  }
}
