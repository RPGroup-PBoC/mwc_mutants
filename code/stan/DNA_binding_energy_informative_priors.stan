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
* normal distributions with user-provided location and 
* scale parameters
*/

#include functions.stan
data {
    // Dimensional parameters
    int<lower=1> N; // Total number of data points.
    
    // Architecdtural parameters
    real R_mu; // mean repressor copy number to be used as the mode
    real R_sig; // standard deviation of repressor copy number
    real<lower=0> Nns; // Number of nonspecific binding sites. 
    
    // Allosteric parameters
    real ep_ai; // Allosteric energy difference in kBT.
    int<lower=1> n_sites; // Number of allosteric binding sites. 
    real<lower=0> Ka_mu; // Inducer dissociation constant to active repressor in units of µM
    real<lower=0> Ka_sig; // Inducer dissociation constant to active repressor in units of µM
    real<lower=0> Ki_mu; // Inducer dissociation constant to inactive repressor in units of µM
    real<lower=0> Ki_sig; // Inducer dissociation constant to inactive repressor in units of µM
    vector<lower=0>[N] c; // Inducer concentration in units of µM
    
    // Observed parameters
    vector[N] fc; // Observed fold-change in gene expression
}

parameters {
   real ep_RA;  // DNA binding energy in units of kBT.
   real Ka; 
   real Ki;
   real R;
   real<lower=0> sigma; // Homoscedastic error

}

model {    
    vector[N] mu;
    
    // Define the weakly informative priors
    ep_RA ~ normal(-12, 6);
    sigma ~ normal(0, .1); 
    
    // Define the informative priors
    R ~ normal(R_mu, R_sig);
    Ka ~ normal(Ka_mu, Ka_sig);
    Ki ~ normal(Ki_mu, Ki_sig);
     
    // Evaluate the likelihood
    for (i in 1:N) {
        mu[i] = fold_change(R, Nns, ep_RA, c[i], log(Ka), log(Ki),
                            ep_ai, n_sites);
        fc[i] ~ normal(mu[i], sigma);
    }
}
