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
    real ref_bohr;
    vector[N] foldchange; // Observed fold-change in gene expression.
}

 
parameters {
    real<lower=0, upper=1> fc_mu; 
    real log_fc_sigma;
}

transformed parameters {
    real fc_sigma = exp(log_fc_sigma);
}
 
model {
    // Define the prior distributions
    fc_mu ~ uniform(0, 1);
    log_fc_sigma ~ normal(0, 10);
    
    // Evaluate the likelihood
    foldchange ~ normal(fc_mu, fc_sigma);
}

generated quantities {
    // Compute the empirical Bohr parameter
    real empirical_bohr = -log((1/fc_mu) - 1); 
    real delta_bohr = ref_bohr - empirical_bohr;
    }