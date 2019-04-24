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
    int<lower=1> N; // Number of measurements
    
    // Architectural parameters 
    real<lower=0> R;  
    real ep_RA;
    real<lower=0> Nns; // Number of non-specific binding sites
    
    // Allosteric parameters
    real<lower=0> Ka;
    real<lower=0> Ki;
    real<lower=0> c[N];
    
    // Observed parameters
    real fc[N]; // Fold-change in gene expression
}

transformed data {
    real ep_a = log(Ka);
    real ep_i = log(Ki);
}

parameters { 
    real ep_AI; // DNA binding energy in kBT
    real<lower=0> sigma; // Homoscedastic error
}

model {
    vector[N] mu;
    
    // Define the priors
    ep_AI ~ normal(0, 5);
    sigma ~ normal(0, 1);
    
    // Evaluate the likelihood
    for (i in 1:N) {
        mu[i] = fold_change(R, Nns, ep_RA, c[i], ep_a, ep_a, ep_AI, 2);
        fc[i] ~ normal(mu[i], sigma);
    } 
    
}

