/* 
* Estimation of DNA Binding Energy from Delta F
* ------------------------------------------------------------------------------
* Author: Griffin Chure
* License: MIT
*
* Description
* ------------------------------------------------------------------------------
* This model infers the DNA binding energy of a given mutant at a given repressor
* copy number by estimating the mean delta F, and computing the mutant DNA 
* binding energy given a reference Bohr parameter. As Delta F itself is an
* inferred quantity in itself,  
*/

data {
    int<lower=1> J; // Unique number of concentrations
    int<lower=1> N; // Number of measurements
    int<lower=1, upper=J> idx[N]; // Identification vector
    vector[N] ref_bohr; // Reference bohr to compute delta bohr
    vector[N] ref_foldchange; // Reference fold-change for calculation of correction factor
    vector[N] foldchange; // Observed fold-change in gene expression.
}
 
parameters {
    real delta_F_mu; // Mean delta F of all measurements
    real<lower=0> sigma; // Standard deviation of delta F for all measurements

    real<lower=0, upper=1> fc_mu[J]; 
    real log_fc_sigma[J]; 
}

transformed parameters {
    real fc_sigma[J] = exp(log_fc_sigma);
    real correction_factor[J];
    real empirical_F[J];
    real delta_F[J];    
    for (i in 1:J) {
    if (fc_mu[i] < fc_sigma[i]) {
        correction_factor[i] = -1 * (ref_bohr[i] + log((1/(ref_foldchange[i] + fc_sigma[i])) - 1));
    }
    else if (1 - fc_mu[i] < fc_sigma[i]) {
        correction_factor[i] = ref_bohr[i] + log((1/(ref_foldchange[i] - fc_sigma[i])) - 1); 
    }
    else {
        correction_factor[i] = 0;
    } 
    empirical_F[i] = -log((1/fc_mu[i]) - 1);
    delta_F[i] = ref_bohr[i] - empirical_F[i] + correction_factor[i];
    }
}
 
model {
    // Define the prior distributions
    delta_F_mu ~ normal(0, 1);
    sigma ~ normal(0, 1); 
    fc_mu ~ uniform(0, 1);
    log_fc_sigma ~ normal(0, 1);
    
    // Evaluate the likelihood
    foldchange ~ normal(fc_mu[idx], fc_sigma[idx]);
    
    for (i in 1:J) {
        delta_F[i] ~ normal(delta_F_mu, sigma); 
    }
}