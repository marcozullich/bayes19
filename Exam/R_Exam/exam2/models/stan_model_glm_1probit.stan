data{
  int<lower=0> N;
  vector[N] x;
  int<lower=0,upper=1> y[N];
}
parameters{
  real alpha;
  real beta;
}
model{
  alpha ~ normal(0,10);
  beta ~ normal(0,2.5);
  y ~ bernoulli(Phi(alpha+beta*x));
}
generated quantities{
  real y_rep[N];
  real log_lik[N];
  
  for(n in 1:N){
    real eta = alpha+beta*x[n];
    y_rep[n] = bernoulli_rng(Phi(eta));
    log_lik[n] = bernoulli_lpmf(y[n] | Phi(eta));
  }
}
