data{
  int n_row;
  vector[n_row] ssb;
  vector[n_row] recruits;
}
parameters{
  real<lower=0> alpha;
  real<lower=0> K;
  real<lower=0> sigma;
  real<lower=0> delta;
}
transformed parameters{
  vector[n_row] mean_recruits;
  
  for(i in 1:n_row){
    mean_recruits[i] = (alpha * ssb[i]^delta)/(1+(ssb[i]^delta)/K);
  }
}
model{
  recruits ~ lognormal(log(mean_recruits)- sigma^2/2, sigma);
  alpha ~ normal(2,2);
  K ~ normal(5000, 15000);
  sigma ~ gamma(2,1);
  delta ~ normal(1,1);
}
generated quantities{
  vector[n_row] log_lik;
  vector[n_row] y_ppd;
  
  for(i in 1:n_row){
    y_ppd[i] = lognormal_rng(log(mean_recruits[i])- sigma^2/2, sigma);
    log_lik[i] = lognormal_lpdf(recruits[i] | log(mean_recruits[i])- sigma^2/2, sigma);
  }
}