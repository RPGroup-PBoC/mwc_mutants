/* 
* Empirical Inference of the Bohr Parameter
* ------------------------------------------------------------------------------
* Author: Griffin Chure
* License: MIT
*
* Description
* ------------------------------------------------------------------------------
* This model infers the fold-change in gene expression for a given set of
* measurements at a given concentration of inducer c. This fold-change is
* restricted to the domain of [0, 1] as the Bohr parameter F is undefined
* outside of these bounds. This model assumes that fold-change is normally 
* distributed with a mean mu and standard deviation sigma. The bound on mu 
* are imposed as 0 and 1 with a uniform prior density in between. This model 
* also calculates the empirical Bohr parameter F from the posterior samples of 
* the mean fold-change.
*/

data {
    int<lower=1> N; // Number of measurements
    int<lower=1> J; // Unique number of measurement sets
    int<lower=1, upper=J> idx[N]; // Identification vector for each measurement
    vector[J] bohr_ref; // Reference Bohr parameter
    vector[N] foldchange; // Observed fold-change in gene expression.
}

parameters {
    vector<upper=0>[J] log_fc_mu;
    vector[J] log_fc_sigma;
}

transformed parameters {
    vector<lower=0, upper=1>[J] fc_mu =exp(log_fc_mu);
    vector<lower=0>[J] fc_sigma = exp(log_fc_sigma);
}

model {
    // Define the prior distributions
    log_fc_mu ~ normal(0, 1); 
    log_fc_sigma ~ normal(-1, 1);
    
    // Evaluate the likelihood
    foldchange ~ normal(fc_mu[idx], fc_sigma[idx]);
}

generated quantities {
    // Compute the empirical Bohr parameter
    vector[J] empirical_bohr;    
    
    // Compute posterior predictive checks to ensure model assumptions are valid. 
    real y_rep[N] = normal_rng(fc_mu[idx], fc_sigma[idx]);
    
    // Evaluate the empirical bohr and delta F
    for (i in 1:J) {
        empirical_bohr[i] = log((1 / fc_mu[i]) - 1); 
    }
}