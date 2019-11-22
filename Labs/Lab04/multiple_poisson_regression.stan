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
  vector<lower=0,upper=1>[N] x2;
  int<lower=0> y[N];
  vector<lower=0>[N] offset;
}

parameters{
  real alpha;
  real beta1;
  real beta2;
  
}

model{
  alpha ~ normal(log(4),1);
  beta1 ~ normal(-.25,1);
  beta2 ~ normal(-.5,1);
  
  y ~ poisson_log(alpha + beta1 * x + beta2 * x2 + offset);
}

generated quantities{
  vector[N] y_rep;
  for (n in 1:N) {
    real eta_n = alpha + beta1 * x[n] + beta2 * x2[n] + offset[n];
    y_rep[n] = poisson_log_safe_rng(eta_n);
  }
}
