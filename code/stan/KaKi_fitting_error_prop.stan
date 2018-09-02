/* 
* Inducer dissociation constant fitting with DNA binding error
* -------------------------------------------------------------
* Author: Griffin Chure
* Date: September 1st, 2018
* License: MIT
* 
* Description
* ---------------------------------------------------------
* Samples the posterior distribution for in the inducer
* dissociation constant to the active (Ka) and inactive (Ki)
* repressor for J unique inducer binding domain mutants. 
* Error in the DNA binding energy is propagated as an informative
* prior. 
*/

data {
    // Dimensional parameters
    int<lower=1> J; // Unique number of inducer binding mutants
    int<lower=1> N; // total number of fold-chnage measurements
    int<lower=1, upper=J> idx[N]; // ID vector for mutants

    // Architectural parameters
    int<lower=0> R_mu; // Mean number of repressors per cell
    int<lower=0> R_sig; // Uncertainty in repressor copy number
    real ep_RA_mu; // Mean value for DNA binding energy in kBT
    real ep_RA_sig; // Uncertainty in DNA binding energy in kBT
    real<lower=0> Nns; // Number of nonspecific binding sites

    // Allosteric parameters
    int n_sites; // Number of alloseteric binding sites
    real ep_AI; // Energy difference betweeen allosteric states in kBT
    vector<lower=0>[N] c; // Effector concentration in µM

    // Observed data
    vector[N] foldchange; // Observed fold-change in gene expression
}


parameters {
    vector<lower=0>[J] Ka;
    vector<lower=0>[J] Ki;
    real<lower=0> R;
    real ep_RA;
    vector<lower=0>[J] sigma;
}

transformed parameters {
    // Perform log transformation of Ka and Ki for better sampling
    vector<lower=-50, upper=50>[J] ep_a; // Inducer binding energy to active repressor in kBT
    vector<lower=-50, upper=50>[J] ep_i; // Inducer binding energy to inactive repreessor in kBT
    ep_a = log(Ka);
    ep_i = log(Ki);
}

model {
    // Instantiate the vector for theoretical prediction
    vector[N] theo;

    // Assign priors
    R ~ normal(R_mu, R_sig);
    ep_RA ~ normal(ep_RA_mu, ep_RA_sig);
    ep_a ~ uniform(-10, 10);
    ep_i ~ uniform(-10, 10);
    sigma ~ normal(0, 1);


    // Evaluate the likelihood
    for (i in 1:N) {
        theo[i] = fold_change(R, Nns, ep_RA, c[i], ep_a[idx[i]], ep_i[idx[i]], ep_AI, n_sites);
        foldchange[i] ~ normal(theo[i], sigma[idx[i]]);
    }
}