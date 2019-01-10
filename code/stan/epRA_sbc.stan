/* 
* Prior Predicitve Checks for DNA Binding Mutants
* -----------------------------------------------
* Author: Griffin Chure
* License: MIT
* 
* Description
* 
#include functions.stan
data {
    // Define the parameters for the simulation based on actual data
    int<lower=0> N; // Number of samples 
    real<lower=0> c[N]; // IPTG concentrations in uM
    real Nns; // Number of nonspecific binding sites
    int<lower=0> n_sites; // Number of allosteric binding sites
    real<lower=0> Ka; // Inducer dissociation constant to active repressor
    real<lower=0> Ki; // Inducer dissociation constant to inactive repressor
    real<lower=0> R;
    real ep_ai; // Energy difference between active and inactive states
    real fc[N]; // Generated ground truth data 
}

transformed data {
    // Log transform the dissocation constants for better sampling
    real ep_a = log(Ka);
    real ep_i = log(Ki);
}

parameters {
    real ep_RA;
    real<lower=0> sigma;
}

model { 
    real mu_sbc[N];
    ep_RA ~ normal(0, 10);
    sigma ~ normal(0, 1);
    for (i in 1:N) {
        mu_sbc[i] = fold_change(R, Nns, ep_RA, c[i],
                               ep_a, ep_i, ep_ai, n_sites);
        fc[i] ~ normal(mu_sbc[i], sigma);
    }
}
