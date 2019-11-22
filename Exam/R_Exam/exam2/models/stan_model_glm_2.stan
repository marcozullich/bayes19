data{
  int<lower=0> N;
  int<lower=0,upper=N> G; //groups
  
  vector[N] x;
  int<lower=0,upper=1> y[N];
  int<lower=1,upper=G> group_mapping[N];
}
parameters{
  vector[G] alpha;
  real beta;
}
model{
  alpha ~ normal(0,10);
  beta ~ normal(0,2.5);
  y ~ bernoulli_logit(alpha[group_mapping]+beta*x);
}
generated quantities{
  real y_rep[N];
  real log_lik[N];
  
  for(n in 1:N){
    real eta = alpha[group_mapping[n]]+beta*x[n];
    y_rep[n] = bernoulli_logit_rng(eta);
    log_lik[n] = bernoulli_lpmf(y[n] | inv_logit(eta));
  }
}
