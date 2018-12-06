/* 
* Double Mutant Parameter Estimation 
* --------------------------------------------------------------
* Author: Griffin Chure
* License: MIT
*
* Description
* -------------------------------------------------------------
* This model samples the posterior probability distributions 
* for Ka, Ki, and ep_RA for all double mutants individually. 
* Due to issues with degeneracy, the allosteric energetic
* difference ep_AI is taken as a delta function. All other 
* parameters (R, Nns, etc) are taken as delta functions at
* their literature values. 
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
    int n_sites; // Number of allosteric inducer binding sites
    real<lower=0> c[N]; // Inducer concentration in µM
     real ep_RA[J]; // DNA binding energy in kBT   
    // Observed parameters
    real fc[N]; // Fold-change in gene expression
}

parameters { 

    real ep_AI[J];
    real<lower=0, upper=1E4> Ka[J]; // Inducer dissociation constant to active repressor
    real<lower=0, upper=1E4> Ki[J]; // Inducer dissociation constant to inactive repressor
    real<lower=0> sigma[J]; // Homoscedastic error
}

transformed parameters {
    // Log transform of Ka/Ki for more efficient sampling. 
    real ep_a[J];
    real ep_i[J];
    for (i in 1:J) {
        ep_a[i] = log(Ka[i]);
        ep_i[i] = log(Ki[i]);
    }
}

model {
    vector[N] mu;
    
    // Define the priors
//    ep_RA ~ normal(0, 10);
    ep_AI ~ normal(0, 10);
    ep_a ~ normal(0, 10);
    ep_i ~ normal(0, 10);
    sigma ~ normal(0, 1);
    
    // Evaluate the likelihood
    for (i in 1:N) {
        mu[i] = fold_change(R, Nns, ep_RA[idx[i]], c[i], ep_a[idx[i]], ep_i[idx[i]],
                          ep_AI[idx[i]], n_sites);
        fc[i] ~ normal(mu[i], sigma[idx[i]]);
    } 
}