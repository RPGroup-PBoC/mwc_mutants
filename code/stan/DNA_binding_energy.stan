/*
* DNA Binding Energy Estimation
* -----------------------------------------------------
* Author: Griffin Chure
* License: CC-BY 4.0
*
* Description
* -----------
* This model saamples the posterior distribution of the DNA binding 
* energy for a set of J unique DNA binding domain mutants of an 
* inducible repressor. All other parameter values (i.e. allosteric energy
* difference, repressor copy number) are taken as delta functions at 
* their literature values. 
*/

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
  vector[N] fc;
}

parameters {
  real ep_RA[J]; // DNA binding energy in units of kBT
  real<lower=0> sigma[J]; //  Homoscedastic error
}

model {
  vector[N] mu;

  // Define the priors. 
  sigma ~ normal(0, 1);
  ep_RA ~ normal(0, 10);

  // Evaluate the likelihood.
  for (i in 1:N) {
    mu[i] = fold_change(R[i], Nns, ep_RA[idx[i]], 0,  0, 0,  ep_ai, n_sites);
    fc[i] ~ normal(mu[i], sigma[idx[i]]);
  }
}