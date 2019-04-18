/*
* Ka/Ki Fitting Model with informative ep_AI
* -----------------------------------------------------
* Author: Griffin Chure
* License: CC-BY 4.0
*
* Description
* ------------------------------------------------------------------------
* This model saamples the posterior distribution of the inducer binding 
* constants to the active and inactive repressor for a set of J unique
* inducer binding domain mutants. This model assumes that each mutant
* has a unique value for the allosteric energy difference (ep_AI),
* inferred from another model in which the value is determined considering
* only the leakiness. All other parameter values (i.e. DNA
* binding energy, allosteric energy difference) are taken as delta
* functions at the literature value. 
*/
#include functions.stan
data { 
  // Dimensional parameters
  int<lower=1> J; // Number of unique mutants.
  int<lower=1> N;  
  int<lower=1, upper=J> idx[N]; // Vector of mutant identifiers
    
  // Architectural parameters
  real<lower=0> R[N]; // Number of repressors
  real<lower=0> Nns; // Number of nonspecific binding sites

  // Allosteric parameters 
  vector<lower=0>[N] c; // Effector concentration.
  real ep_RA; // Binding energy in kBT. 
  int<lower=1> n_sites; // Number of allosteric sites.  

  // Observed parameters.
  vector[N] fc;
  }

parameters {
  real<lower=0> Ka[J]; // log transform of Ka
  real<lower=0> Ki[J]; // log transform of Ki
  real ep_AI[J];
  real<lower=0> sigma[J]; //  Homoscedastic error
}

transformed parameters {
  real ep_a[J] = log(Ka); 
  real ep_i[J] = log(Ki);
}

model {
  vector[N] mu;
    
  // Define the priors. 
  sigma ~ normal(0, 0.1);
  Ka ~ lognormal(2, 2);
  Ki ~ lognormal(0, 2);
  ep_AI ~ normal(0, 5);

  for (i in 1:N) {
    mu[i] = fold_change(R[i], Nns, ep_RA, c[i], ep_a[idx[i]], ep_i[idx[i]], ep_AI[idx[i]], n_sites);
    fc[i] ~ normal(mu[i], sigma[idx[i]]);
  }
}