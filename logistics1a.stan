data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real<lower=0, upper=1> alpha;
  real<lower=0, upper=1>  beta;
  real<lower=0, upper=1> sigma;
  real<lower=1, upper=6> delta;
  real<lower=0> eps;
}
model {
  for (n in 1:N)
    y[n] ~ normal((alpha * pow(x[n],sigma))/(delta + beta * pow(x[n],sigma)), eps); 
    increment_log_prob(-log(eps));     //log prior for p(eps) proportion to 1/eps
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N){
    log_lik[n] <- normal_log(y[n], (alpha * pow(x[n],sigma))/(delta + beta * pow(x[n],sigma)), eps);
  }
}