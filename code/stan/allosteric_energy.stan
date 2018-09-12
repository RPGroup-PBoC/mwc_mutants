/*
* Inference of ep_AI
* --------------------------------------------------------------
*  Author: Griffin Chure
*  License: MIT
* 
* Description
* --------------------------------------------------------------
* This model samples the posterior distribuiton for the energetic
* difference between the two allosteric states given a set of 
* leakiness measurements at a single repressor copy number. 
* The results from this model can then be used to inform the
* inference of the inducer dissociation constants Ka/Ki for the 
* full titration 
*/ 
data {
    // Dimensional parameters
    int<lower=1> J; // Number of unique inducer binding mutants
    int<lower=1> N; // Total number of data points. 
    int<lower=1, upper=J> idx[N]; // Identification vector for mutants

    // Architectural parameters
    vector<lower=0>[N] R; // Average number of repressors per cell
    real ep_RA; // WT DNA Binding energy in kBT.
    real<lower=0> Nns; // Number of nonspecific binding sites. 

    // Measured parameters
    vector<lower=-0.2, upper=1.2>[N] fc; 
}

parameters {
   real ep_ai[J]; // Allosteric energy difference in kBT
   real<lower=0> sigma[J]; // Homoscedastic error
}

model {
    vector[N] mu;

    // Priors
    ep_ai ~ normal(0, 10);
    sigma ~ normal(0, 1);

    // Compute the expectation value
    for (i in 1:N) {
        mu[i] = fold_change(R[i], Nns, ep_RA, 0, 1, 1, ep_ai[idx[i]], 1);       
        fc[i] ~ normal(mu[i], sigma[idx[i]]);
    }
}

    
