/* 
* Allosteric Energy Difference Fitting
* --------------------------------------------------------------
* Author: Griffin Chure
* License: MIT
*
* Description
* -------------------------------------------------------------
* This model samples the posterior probability distributions 
* for the energetic difference between allosteric states. 
* All other parameters (R, Nns, etc) are taken as delta functions at
* their provided values. Only the leakiness is considered for this
* fitting.
*/
#include functions.stan

data {
    // Dimensional parameters
    int<lower=1> J; // Number of unique mutants
    int<lower=1> N; // Number of measurements
    int<lower=1, upper=J> idx[N]; // Identification vector
    
    // Architectural parameters 
    real<lower=0> Nns; // Number of non-specific binding sites
    
    // Allosteric parameters
    real<lower=0> R;  
    real ep_RA[J];
    real ep_AI[J];
    
    // Observed parameters
    real fc[N]; // Fold-change in gene expression
}

parameters { 
    real ep_int[J]; // DNA binding energy in kBT
    real<lower=0> sigma[J]; // Homoscedastic error
}

model {
    vector[N] mu;
    
    // Define the priors
    ep_int ~ normal(0, 10);
    sigma ~ normal(0, 1);
    
    // Evaluate the likelihood
    for (i in 1:N) {
        mu[i] = fold_change(R, Nns, ep_RA[idx[i]], 0, 1, 1,
                          ep_AI[idx[i]] + ep_int[idx[i]], 2);
        fc[i] ~ normal(mu[i], sigma[idx[i]]);
    } 
}