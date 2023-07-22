data {
  int N; // number of datapoints
  vector[N] y; // stream density
  int Nstreams;
  int<lower=0,upper=Nstreams> stream_id[N];
}

parameters {
  real beta;
  real<lower=0> sigma;
  vector[Nstreams] bi; // stream-effects
  real<lower=0> sigma_stream; // variance in stream effects
}

transformed parameters {
  vector[N] mu;
  for(i in 1:N)
  {
    mu[i] = beta + bi[stream_id[i]];
  }
}

model {
  beta ~ normal(50,100);
  bi ~ normal(0,sigma_stream);
  sigma ~ gamma(2,1);
  sigma_stream ~ gamma(2,1);
  y ~ normal(mu,sigma);
}

generated quantities {
  vector[N] log_lik;
  vector[N] y_ppd;
  vector[Nstreams] stream_effect;
  for(i in 1:Nstreams)
  {
    stream_effect[i] = beta + bi[i];
  }
  for(i in 1:N)
  {
    y_ppd[i] = normal_rng(mu[i],sigma);
    log_lik[i] = normal_lpdf(y[i] | mu[i], sigma);
  }
}

