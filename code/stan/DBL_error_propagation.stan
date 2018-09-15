/*
* Double Mutant Error Propagation
* --------------------------------------------------------------------
* Author: Griffin Chure
* License: MIT
*
* Description
* --------------------------------------------------------------------
* This model samples the posterior distribution of DNA binding
* energy as well as the inducer dissociation constants for 
* double mutants of an inducible repressor. This model applies
* highly informative priors for the various fitted parameter values 
* from the single mutants with uniform distributions between the 
* lower and upper bound of the 95% credible region for each parameter. 
* The purpose of this restrictive sampling is to give a sense of the 
* uncertainty in the predictions generated using the values from the 
* single mutants. 
*/

data {
    // Dimensional parameters. 
    int<lower=1> J_DNA; // Unique number of DNA binding domain mutants
    int<lower=1> J_IND; // Unique number of inducer binding domain mutants
    int<lower=1> J; // Unique number of double mutants
    int<lower=1> N; // Total number of measurements
    int<lower=1, upper=J_DNA> DNA_idx[N]; // Vector of DNA mutant identifiers
    int<lower=1, upper=J_IND> IND_idx[N]; // Vector of inducer mutant identifiers
    int<lower=1, upper=J> idx[N]; // Vector of double mutant identifiers

    // Architectural parameters. 
    real<lower=0> Nns; // Number of nonspecific binding sites. 
    real<lower=0> R; // Average number of repressors per cell.

    // Allosteric parameters
    int<lower=0> n_sites; // Number of allosteric binding sites.
    real<lower=0> c[N]; // Inducer concentration in µM.
    real ep_AI_mu[J_IND]; // Center of informative prior for allosteric energy difference in kBT.
    real ep_AI_sig[J_IND]; // Variance of informative for allosteric energy difference in kBT.

    // Informative prior bounds
    real ep_RA_mu[J_DNA]; // Center of informative prior for ep_RA
    real<lower=0> ep_RA_sig[J_DNA]; // Variance for informative prior
    real<lower=0> Ka_mu[J_IND]; // in µM
    real<lower=0> Ka_sig[J_IND]; //  in µM 
    real<lower=0> Ki_mu[J_IND]; // in µM 
    real<lower=0> Ki_sig[J_IND]; // in µM

    // Observed data
    vector<lower=-0.2, upper=1.2>[N] fc; // Observed fold-change in gene expression
}

parameters {
    real<upper=0> ep_RA[J_DNA];
    real<lower=0> Ka[J_IND];
    real<lower=0> Ki[J_IND];
    real<lower=0> ep_AI[J_IND];
    real<lower=0> sigma[J]; // Homoscedastic error
}

transformed parameters {
     // Log transform of inducer dissociation constants for more efficient sampling. 
     real ep_a[J_IND];
     real ep_i[J_IND];
     for (i in 1:J_IND) {
         ep_a[i] = log(Ka[i]);
         ep_i[i] = log(Ki[i]);
     }
 }

model {
    // Instantiate a vector of the theoretical value. 
    vector[N] mu;

    // Define the priors as informative uniform. 
    for (i in 1:J_DNA) {
        ep_RA[i] ~ normal(ep_RA[i],  ep_RA_sig[i]);
    }

    for (i in 1:J_IND) {
        Ka[i] ~ normal(Ka_mu[i], Ka_sig[i]);
        Ki[i] ~ normal(Ki_mu[i], Ki_sig[i]); 
        ep_AI[i] ~ normal(ep_AI_mu[i], ep_AI_sig[i]); 
    }

    // Define prior for homoscedastic errors. 
    sigma ~ normal(0, 1);

    // Evaluate the likelihood. 
    for (i in 1:N) {
        mu[i] = fold_change(R, Nns, ep_RA[DNA_idx[i]], c[i], ep_a[IND_idx[i]], ep_i[IND_idx[i]],
                           ep_AI[IND_idx[i]],n_sites);
        fc[i] ~ normal(mu[i], sigma[idx[i]]);
    }
}