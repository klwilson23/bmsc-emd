data {
  int N; // the number of datapoints
  int Ncovar; // the number of coefficients to estimate
  int lice[N]; // the count of lice on each burb
  matrix[N,Ncovar] Xvar; // the coefficients
  int N_aviary; // the number of aviaries - one of our random effects
  int N_burbs; // the number of breeds
  int<lower=0> aviary_id[N]; // the group indicat for each aviary
  int<lower=0> burb_id[N]; // the indicator for each breed
  int burb_ave[N_aviary]; // the indicator for which breed is in each aviary
  int N_predict; // the number of rows of a prediction matrix
  matrix[N_predict,Ncovar] x_predict; // a matrix of predictions for easy plotting
}

parameters {
  vector[Ncovar] betas;
  real<lower=0> sigma_aviary;
  real<lower=0> sigma_burb;
  vector[N_burbs] burb_eff;
  vector[N_aviary] aviary_eff;
  real<lower=0> phi;
}

transformed parameters {
  vector[N] log_mu;
  vector[N] mu;
  for(i in 1:N)
  {
      log_mu[i] = (Xvar[i]*betas);
      mu[i] = exp(log_mu[i] + aviary_eff[aviary_id[i]]);
  }
}

model {
  betas ~ normal(0,2);
  lice ~ neg_binomial_2(mu,phi);
  burb_eff ~ normal(0,sigma_burb);
  phi ~ inv_gamma(0.1,0.1);
  sigma_burb ~ gamma(2,1);
  sigma_aviary ~ gamma(2,1);
  
  for(i in 1:N_aviary)
  {
    aviary_eff[i] ~ normal(burb_eff[burb_ave[i]],sigma_aviary);
  }
}

generated quantities {
  vector[N] y_ppd;
  vector[N] log_lik;
  vector[N_predict] new_mu;
  for(i in 1:N)
  {
    y_ppd[i] = neg_binomial_2_rng(mu[i],phi);
    log_lik[i] = neg_binomial_2_lpmf(lice[i]|mu[i],phi);
  }
  for(i in 1:N_predict)
  {
    new_mu[i] = exp(x_predict[i]*betas);
  }
}
