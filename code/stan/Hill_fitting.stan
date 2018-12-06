/* 
* Hill Function Inference 
* ------------------------------------------
* Author: Griffin Chure 
* License: MIT
* 
* Description
* -----------------------------------------
* This model samples the posterior probability distribution 
* for the various parameters of a Hill function to a 
* single titration curve. 
* 
* The generalized equation for a Hill function is
* fold-change = a + b * ((c / K)^n / 1 + (c / K)^n)
* where a is the leakiness, b is the dynamic range,
* c is the inducer concentration, K is the effective 
* inducer repressor dissociation constant, and n
* is the HIll coefficient. 
*/

data  {
    
    // Dimensional parameters
    int<lower=1> J; // Unique number of mutants
    int<lower=1> N; // Total number of measurements. 
    int<lower=1, upper=J> idx[N]; // Identification vector
    
    // Observed parameters
    real<lower=0> c[N]; // Inducer concentration
    real foldchange[N]; // Fold-change in gene expression
}

parameters { 
    real<lower=0, upper=1> a[J];
    real<lower=0, upper=1> b[J];
    real<lower=0> K[J];
    real<lower=0> n[J]; 
    real<lower=0> sigma[J];
 
}
 
model {
    // Define a vector for the theoretical fold-change
    vector[N] mu;
    
    // Define the priors. Assume a half normal for zero bounded and normal for all others.
    a ~ normal(0, 1);
    b ~ normal(0, 1);
    K ~ lognormal(0, 1);
    n ~ normal(0, 1);
    sigma ~ normal(0, 1);
    
    // Evaluate the likelihood 
    for (i in 1:N) {
        mu[i] = a[idx[i]] + b[idx[i]] * ((c[i] / K[idx[i]])^n[idx[i]] / (1 + (c[i] / K[idx[i]])^n[idx[i]]));
        foldchange[i] ~ normal(mu[i], sigma[idx[i]]);
    }
}