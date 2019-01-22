/* 
* Parameter Estimation via Delta Delta F
* --------------------------------------
* Author: Griffin Chure
* License: MIT
* 
* Description 
* -------------------------------------
* This model performs parameter estimation of 
* Ka, Ki, epAI, and epRA from the measurements
* given the empirical measurement of delta delta F.
*/

#include functions.stan

data {
    int<lower=1> N; // Total number of measurements 

    // Wild-type parameters. 
    real<lower=0> wt_Ka; // In uM
    real<lower=0> wt_Ki; // In uM
    int<lower=1> n_sites; // Number of allosteric binding sites. 
    vector<lower=0>[N] c; // IPTG concentration in uM
    real wt_ep_AI; // Energy difference between allosteric states in kT
    real wt_ep_RA; // DNA binding energy in kT

    // Observed parameters
    vector[N] ddF; // Difference between empirical  
}

transformed data {
    // Perform log transforms of wild type inducer dissociation constants. 
    real wt_ep_a = log(wt_Ka);
    real wt_ep_i = log(wt_Ki);

    // Compute the wild-type P_act. 
    vector<lower=0, upper=1>[N] wt_pact;
    for (i in 1:N) {
        wt_pact[i] = prob_act(c[i], wt_ep_a, wt_ep_i, wt_ep_AI, n_sites);
    }
}

parameters {
    // Define all biophysical parameters of the mutant. 
    real ep_RA; // DNA binding energy of mutant
    real ep_AI; // Allosteric energy difference. 
    real<lower=0> Ka; // in uM
    real<lower=0> Ki; // in uM
    real<lower=0, upper=1> sigma;
}

transformed parameters {
    // Perform log transform of estimated dissociation constants.
    real ep_a = log(Ka);
    real ep_i= log(Ki);

    // Compute the theoretical ddf. 
    vector<lower=0, upper=1>[N] pact;
    vector[N] theo_ddF;
    real ep_RA_diff = ep_RA - wt_ep_RA;

    for (i in 1:N) {
        pact[i] = prob_act(c[i], ep_a, ep_i, ep_AI, n_sites);
        theo_ddF[i] = log(wt_pact[i] / pact[i]) - ep_RA_diff;
    }
}

model {
    // Define the priors.
    ep_RA ~ normal(-12, 6);
    ep_AI ~ normal(0, 5);
    ep_a ~ normal(0, 5);
    ep_i ~ normal(0, 5);
    sigma ~ normal(0, 0.1);

    // Compute the likelihood    
    ddF ~ normal(theo_ddF, sigma);
}