/* 
* DNA Binding Energy Estimation from Leakiness Measurement
* --------------------------------------------------------
* Author: Griffin Chure
* Date: August 8, 2018
* License: MIT
* 
* Description
* ---------------------------------------------------------
* Samples the posterior distribution for the DNA binding 
* energy for J unique DNA binding domain mutants with proper
* propagation of error from repressor copy number determination.
* In this model, only the measurements of the leakiness are 
* considered for estimation. 
*/

data {
    // Dimensional parameters
    int<lower=1> J; // Unique number of DNA mutants. 
    int<lower=1> J_rep; // Total number of unique repressor copy numbers
    int<lower=1> N; // Total number of fold-change measurements 
    int<lower=1, upper=J> idx[N]; // ID vector for mutants
    int<lower=1, upper=J_rep> rep_idx[N]; // ID vector for repressor copy number


    // Architectural parameters
    vector<lower=0>[J_rep] R_mu; // Mean repressor copy number
    vector<lower=0>[J_rep] R_sig; // Uncertainty in repressor copy number
    real Nns; // Number of nonspecific binding sites

    // Allosteric parameters
    int n_sites; // Number of allosteric binding sites
    real ep_AI; // Energy difference between allosteric states in kBT

    // Observed variables.
    vector[N] foldchange; // Observed fold-change in gene expression
}

parameters {
    vector<lower=-30, upper=30>[J] ep_RA;  // DNA binding energy for each mutant in kBT
    vector<lower=0>[J_rep] R;
    vector<lower=0>[J] sigma; // Homoscedastic error in fold-change measurement
}

transformed parameters {
    vector[N] foldchange_log;
    foldchange_log = log(foldchange);
}
model {
    // Instantiate the vector for theoretical prediction
    vector[N] theo;

    // Assign priors    
    R ~ normal(R_mu, R_sig);
    ep_RA ~ normal(0, 10);
    sigma ~ normal(0, 1);

    // Compute the likelihood
    for (i in 1:N) {
        theo[i] = fold_change(R[rep_idx[i]], Nns, ep_RA[idx[i]], 0, 1, 1, ep_AI, n_sites);
        foldchange_log[i] ~ normal(log(theo[i]), sigma[idx[i]]);
    }
}

