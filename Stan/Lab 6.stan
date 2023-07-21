data {
  int N;
  vector[N] density;
  int Nstreams;
  int<lower=0,upper=Nstreams> stream_id[N];
  
}

parameters {
  real beta;
  vector[Nstreams] stream_effect;
  real<lower=0> sigma_stream;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] mu;
  for(i in 1:N)
  {
    mu[i] = beta + stream_effect[stream_id[i]];
  }
}

model {
  beta ~ normal(50,100);
  stream_effect ~ normal(0,sigma_stream);
  y ~ normal(mu,sigma);
  sigma_stream ~ gamma(2,1);
  sigma ~ gamma(2,1);
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_ppd;
  
  for(i in 1:N){
    y_ppd[i] = normal_rng(mu[i], sigma);
    log_lik[i] = normal_lpdf(y[i] | mu[i], sigma);
  }
}
