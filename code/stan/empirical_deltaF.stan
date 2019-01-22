/* 
* Inference of the empirical delta F
* ---------------------------------------------------
* Author: Griffin Chure
* License: MIT
*
* Description 
* ----------------------------------------------------
* This Stan program performs inference of the mean
* fold-change at a range of IPTG concentrations given 
* a set of measurements. This value, bounded between 0
* and 1, is then used to compute the measured value of
* the Bohr parameter F. From this, the empirical delta 
* F is calculated given a reference value. 
*/ 

data {
    // Dimensional parameters
    int<lower=1> N; // Total number of measurements. 
    int<lower=1> J; // Number of unique IPTG concentrations
    int<lower=1, upper=J> idx[N]; // Identification vector for each measurement. 
    
    // Architectural parameters. 
    real F0; // Reference Bohr parameter.
    real deltaF_ref[J]; // Delta F of reference strain at each IPTG concentration
     
    // Observed parameters
    vector[N] fc; // Observed value of fold-change.     
}

parameters {
    vector<lower=0, upper=1>[J] fc_mu; // Mean fold-change at each IPTG concentration
    real<lower=0, upper=1> sigma; // Homoscedastic error in fold-change measurement.
}

model {
    // Define the priors
    fc_mu ~ uniform(0, 1);
    sigma ~ normal(0, 0.1);    
    
    // Evaluate the likelihood. 
    fc[idx] ~ normal(fc_mu[idx], sigma);
}

generated quantities {
    // Compute the measured Bohr parameter and empirical delta F.
    vector[J] measured_F;
    vector[J] empirical_delF;
    vector[J] deldelF;
    for (i in 1:J) {
        measured_F[i] = log((1 / fc_mu[i]) - 1);
        empirical_delF[i] = F0 - measured_F[i];
        deldelF[i] = empirical_delF[i] - deltaF_ref[i];
    }
}