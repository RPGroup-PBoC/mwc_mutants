/*
* Ka/Ki Fitting Model
* -----------------------------------------------------
* Author: Griffin Chure
* License: CC-BY 4.0
*
* Description
* -----------
* This model saamples the posterior distribution of the inducer binding 
* constants to the active and inactive repressor for a set of J unique
* inducer binding domain mutants. All other parameter values (i.e. DNA
* binding energy, allosteric energy difference) are taken as delta
* functions at the literature value. 
*/

data { 
  // Dimensional parameters
  int<lower=1> J; // Number of unique mutants.
  int<lower=1> N; // Total number of data points.
  int<lower=1, upper=J> idx[N]; // Vector of mutant identifiers

  // Architectural parameters
  vector<lower=0>[N] R; // Number of repressors
  real<lower=0> Nns; // Number of nonspecific binding sites

  // Allosteric parameters 
  vector<lower=0>[N] c; // Effector concentration.
  real ep_RA; // Binding energy in kBT. 
  real ep_AI; // Allosteric energy difference 
  int<lower=1> n_sites; // Number of allosteric sites.  

  // Observed parameters.
  vector<lower=0, upper=1.2>[N] fc;
  }

parameters {
  real<lower=0, upper=5000> Ka[J]; // Active repressor inducer dissociation constant
  real<lower=0, upper=5000> Ki[J]; // Inactive repressor inducer dissociation constant
  real<lower=0> sigma[J]; //  Homoscedastic error
}

transformed parameters {
  real ep_a[J]; 
  real ep_i[J];
  ep_a = log(Ka);
  ep_i = log(Ki);
}

model {
  vector[N] mu;

  // Define the priors. 
  sigma ~ normal(0, 1);
  ep_a ~ normal(0, 10);
  ep_i ~ normal(0, 10);

  for (i in 1:N) {
    mu[i] = fold_change(R[i], Nns, ep_RA, c[i], ep_a[idx[i]], ep_i[idx[i]], ep_AI, n_sites);
    fc[i] ~ normal(mu[i], sigma[idx[i]]);
  }
}