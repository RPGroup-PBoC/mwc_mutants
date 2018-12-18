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
    real<lower=0> Ka[J];
    real<lower=0> Ki[J];
    real<lower=0> c[N];
    
    // Observed parameters
    real fc[N]; // Fold-change in gene expression
}

transformed data {
    real ep_a[J] = log(Ka);
    real ep_i[J] = log(Ki);
}

parameters { 
    real ep_AI[J]; // DNA binding energy in kBT
    real<lower=0> sigma[J]; // Homoscedastic error
}

model {
    vector[N] mu;
    
    // Define the priors
    ep_AI ~ normal(0, 10);
    sigma ~ normal(0, 1);
    
    // Evaluate the likelihood
    for (i in 1:N) {
        mu[i] = fold_change(R, Nns, ep_RA[idx[i]], c[i], ep_a[idx[i]], ep_a[idx[i]], ep_AI[idx[i]], 2);
        fc[i] ~ normal(mu[i], sigma[idx[i]]);
    } 
}