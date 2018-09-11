/* 
* DNA Binding Energy Estimation from Induction Profile
* --------------------------------------------------------
* Author: Griffin Chure
* License: MIT

* Description
* -----------
* This model samples the posterior distribution of the DNA 
* binding energy for a set of J unique induction profiles. 
* All other parameter values (i.e. allosteric energy 
* difference, Ka, Ki, etc) are taken as
* delta functions at their literature value. 
*/

data {
    // Dimensional parameters
    int<lower=1> J; // Number of uniqueinduction profiles.
    int<lower=1> N; // Total number of data points.
    int<lower=1, upper=J> idx[N]; // Identification vector for induction profile. 
    
    // Architecdtural parameters
    vector<lower=0>[N] R; // Average number of repressors per cell.
    real<lower=0> Nns; // Number of nonspecific binding sites. 
    
    // Allosteric parameters
    real ep_ai; // Allosteric energy difference in kBT.
    int<lower=1> n_sites; // Number of allosteric binding sites. 
    real<lower=0> Ka; // Inducer dissociation constant to active repressor in units of µM
    real<lower=0> Ki; // Inducer dissociation constant to inactive repressor in units of µM
    vector<lower=0>[N] c; // Inducer concentration in units of µM
    // Observed parameters
    vector[N] fc; // Observed fold-change in gene expression
}

parameters {
    real ep_RA[J];  // DNA binding energy in units of kBT.
    real<lower=0> sigma[J]; // Homoscedastic error
}

model {
    vector[N] mu; // Storage vector for theoretical value.
    
    // Define the priors
    ep_RA ~ normal(0, 10);
    sigma ~ normal(0, 1); 
    
    // Evaluate the likelihood
    for (i in 1:N) {
        mu[i] = fold_change(R[i], Nns, ep_RA[idx[i]], c[i], log(Ka), log(Ki),
                            ep_ai, n_sites);
        fc[i] ~ normal(mu[i], sigma[idx[i]]);
    }
}
    
