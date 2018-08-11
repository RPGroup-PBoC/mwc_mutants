/*
* Error Propagation for DNA Mutant Predictions
* --------------------------------------------
* Author: Griffin Chure
* Date: 10 August 2018
* License: MIT
*
* Description
* --------------------------------------------
* This model performs the propagation of uncertainty
* for testing the predictions for fitting the DNA binding
* energy to the single-mutants in the DNA binding domain. 
* This model samples the model using highly informative 
* priors implemented as normal distributions.
*/

data {
    // Dimensional parameters
    int<lower=1> J; // Nunique number of DNA mutants
    int<lower=1> J_rep; // Total number of unique repressor copy numbers
    int<lower=1> N; // Total number of fold-change measurements
    int<lower=1, upper=J> idx[N]; // ID vector for mutants
    int<lower=1, upper=J_rep> rep_idx[N]; // ID vector for repressor copy number
    
    // Architectural parameters
    real<lower=0> R_mu[J_rep]; // Mean repressor copy number
    real<lower=0> R_sig[J_rep]; // Uncertainty in repressor copy number
    vector[J] ep_RA_mu; // Mean DNA binding energy
    vector[J] ep_RA_sig; // Uncertainty in DNA binding energy
    real Nns; // Number of nonspecific binding sties
    
    // Allosteric parameters
    real<lower=0> Ka_mu; // Mean inducer dissociation constant to active repressor in µM
    real<lower=0> Ka_sig; // Uncertainty in active repressor inducer binding constant in µM
    real<lower=0> Ki_mu; // Mean inducer dissociation constant to inactive repressor in µM
    real<lower=0> Ki_sig; // Uncertainty in inactive repressor inducer binding constant in µM
    int n_sites; // Number of allosteric binding sites
    real ep_AI; // Energy difference between allosteric sites in kBT
    vector<lower=0>[N] c; // Inducer concentration in µM
    
    // Observed variables
    vector[N] foldchange; // Fold-change in gene expression
}

parameters {
    vector[J] ep_RA; 
    real<lower=0> Ka; 
    real<lower=0> Ki;
    vector<lower=0>[J_rep] R;
    vector<lower=0>[J] sigma; // Homoscedastic error
}


transformed parameters {
    // Perform log transformation of inducer dissociation constants for better sampling
    real ep_a; // Active repressor inducer binding energy in kBT 
    real ep_i; // Inactive repressor inducer binding energy in kBT
    ep_a = log(Ka);
    ep_i = log(Ki);
}

model {
    // Instantiate the vector for theoretical value
    vector[N] theo;
    
    // Define the informative priors
    ep_RA ~ normal(ep_RA_mu, ep_RA_sig);
    R ~ normal(R_mu, R_sig);
    Ka ~ normal(Ka_mu, Ka_sig);
    Ki ~ normal(Ki_mu, Ki_sig);
    sigma ~ normal(0, 1);
    
    // Iterate through the data and compute the likelihood. 
    for (i in 1:N) {
        theo[i] = fold_change(R[rep_idx[i]], Nns, ep_RA[idx[i]], c[i], ep_a, ep_i, 
                             ep_AI, n_sites);
        foldchange[i] ~ normal(theo[i], sigma[idx[i]]);
    }
}