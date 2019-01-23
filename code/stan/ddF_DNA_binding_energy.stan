/* 
* Parameter Estimation via Delta Delta F
* ------------------------------------------------------------------------------
* Author: Griffin Chure
* License: MIT
* 
* Description 
* ------------------------------------------------------------------------------
* This model performs parameter estimation of Ka, Ki, epAI, and epRA from the 
* measurements given the empirical measurement of delta delta F.
*/
#include functions.stan
data {
    int<lower=1> N; // Total number of measurements 
    real wt_ep_RA; // DNA binding energy in kT
    vector[N] ddF; // Difference between empirical  
}

parameters {
    // Define all biophysical parameters of the mutant. 
    real ep_RA; // DNA binding energy of mutant
    real<lower=0, upper=1> sigma;
}

transformed parameters {
    // Perform log transform of estimated dissociation constants.
    for (i in 1:N) {
        pact[i] = prob_act(c[i], ep_a, ep_i, ep_AI, n_sites);
        theo_ddF[i] = log(wt_pact[i] / pact[i]) - ep_RA_diff;
    }
}

model {
    // Define the priors.
    ep_RA ~ normal(-12, 6);
    sigma ~ normal(0, 0.1);

    // Compute the likelihood    
    ddF  ~ normal(ep_RA - wt_ep_RA, sigma);
}