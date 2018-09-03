data { 
  // Dimensional parameters
  int J; // Number of unique mutants
  int N; // Total number of data points.
  int idx[N]; // Trial idx

  // Architectural parameters
  vector<lower=0>[N] R; // Number of repressors
  real<lower=0> Nns; // Number of nonspecific binding sites

  // Allosteric parameters 
  real ep_ai; // Allosteric energy difference 
  int n_sites; // Number of allosteric binding sites

  // Observed parameters.
  vector<lower=-0, upper=1.25>[N] fc;
}

parameters {
  real ep_RA[J]; // DNA binding energy in units of kBT
  real<lower=0> sigma[J]; //  Homoscedastic error
}

transformed parameters {
  vector[N] log_fc;
  log_fc = log(fc);
}

model {
  vector[N] mu;

  // Define the priors. 
  sigma ~ normal(0, 1);
  ep_RA ~ normal(0, 10);

  for (i in 1:N) {
    mu[i] = fold_change(R[i], Nns, ep_RA[idx[i]], 0,  0, 0,  ep_ai, n_sites);
    log_fc ~ normal(log(mu[i]), sigma[idx[i]]);
  }
}