
data{
  int N;          // number of voters
  int vote[N];    // vote: 0 (Clinton), 1 (Bush)
  int income[N];  // 1-5 income scale
}
parameters{
  real alpha;    // intercept
  real beta;     // income coefficient
}
model{
  //alpha ~ normal(0,10);
  //beta ~ normal(0,2.5);
  //for(i in 1:N){
  //  vote[i] ~ bernoulli(Phi(alpha + beta*income[i]));
  //}
  target += normal_lpdf(alpha | 0, 10);
  target += normal_lpdf(beta | 0, 2.5);
  for(i in 1:N){
    target += bernoulli_lpmf(vote[i]  | Phi(alpha + income[i]*beta));
  }
   
}

generated quantities{
  vector[N] y_rep;
  for (n in 1:N) {
    real eta_n = alpha + beta * income[n];
    y_rep[n] = bernoulli_rng(Phi(eta_n));
  }
}
