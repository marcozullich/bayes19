functions{
  int poisson_log_safe_rng(real eta){
    real pois_rate = exp(eta);
    if(pois_rate >= exp(20.79)){
      return -9;
    }
    return poisson_rng(pois_rate);
  }
}

data {
  int<lower=1> N;
  vector<lower=0>[N] x;
  int<lower=0> y[N];
}

parameters{
  real alpha;
  real beta;
  
}

model{
  alpha ~ normal(log(4),1);
  beta ~ normal(-.25,1);
  
  y ~ poisson_log(alpha + beta * x);
}

generated quantities{
  vector[N] y_rep;
  for (n in 1:N) {
    real eta_n = alpha + beta * x[n];
    y_rep[n] = poisson_log_safe_rng(eta_n);
  }
}
