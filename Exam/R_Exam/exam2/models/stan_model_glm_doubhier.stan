data{
  int<lower=0> N;
  int<lower=0,upper=N> G1; //groups1
  int<lower=0,upper=N> G2; //groups1
  
  vector[N] x;
  int<lower=0,upper=1> y[N];
  int<lower=1,upper=G1> group_mapping1[N];
  int<lower=1,upper=G2> group_mapping2[N];
}
parameters{
  vector[G1] alpha1;
  vector[G2] alpha2;
  real beta;
}
model{
  alpha1 ~ normal(0,1);
  alpha2 ~ normal(0,1);
  beta ~ normal(0,2.5);
  y ~ bernoulli_logit(alpha1[group_mapping1]+alpha2[group_mapping2]+beta*x);
}
generated quantities{
  real y_rep[N];
  real log_lik[N];
  
  for(n in 1:N){
    real eta = alpha1[group_mapping1[n]]+alpha2[group_mapping2[n]]+beta*x[n];
    y_rep[n] = bernoulli_logit_rng(eta);
    log_lik[n] = bernoulli_lpmf(y[n] | inv_logit(eta));
  }
}
