/*
* Double Mutant Error Propagation
* -----------------------------------------------------
* Author: Griffin Chure
* License: CC-BY 4.0
*
* Description
* -----------
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
    real<lower=0> Nns; // Njumber of nonspecific binding sites. 
    real<lower=0> R; // Average number of repressors per cell.

    // Allosteric parameters
    int<lower=0> n_sites; // Number of allosteric binding sites.
    real<lower=0> c; // Inducer concentration in µM.
    real ep_AI; // Allosteric energy difference in kBT.

    // Informative prior bounds
    real<upper=0> ep_RA_lower[J_DNA]; // Lower bound of the DNA binding energy in kBT.
    real<upper=0> ep_RA_upper[J_DNA]; // Upper bound of the DNA binding energy in kBT.
    real<lower=0> Ka_upper[J_IND]; // Upper bound of the active repressor inducer dissociation constant in µM
    real<lower=0> Ka_lower[J_IND]; // Lower bound of the active repressor inducer dissociation constant in µM 
    real<lower=0> Ki_upper[J_IND]; // Upper bound of the inactive repressor inducer dissociation constant in µM
    real<lower=0> Ki_lower[J_IND]; // Lower bound of the inactive repressor inducer dissociation constant in µM 


    // Observed data
    vector<lower=0, upper=1.2> fc[N]; // Observed fold-change in gene expression
}

parameters {
    real<upper=0> ep_RA[J_DNA];
    real<lower=0> Ka[J_IND];
    real<lower=0> Ki[J_IND];
    real<lower=0> sigma[J]; // Homoscedastic error
}

transformed parameters {
    // Log transform of inducer dissociatoin constants for more efficient sampling. 
    real ep_a[J_IND];
    real ep_i[J_IND];
    ep_a = log(Ka);
    ep_i = log(Ki);
}

model {
    // Instantiate a vectof or the theoretical value. 
    vector<lower=0, upper=1> mu[N];

    // Define the priors as informative uniform. 
    for (i in 1:J_DNA) {
        ep_RA[i] ~ uniform(ep_RA_lower[i],  ep_RA_upper[i]);
    }

    for (i in 1:J_IND) {
        Ka[i] ~ uniform(Ka_lower[i], Ka_upper[i]);
        Ki[i] ~ uniform(Ki_lower[i], Ki_upper[i]);
    }

    // Define prior for homoscedastic errors. 
    sigma ~ normal(0, 1);

    // Evaluate the likelihood. 
    for (i in 1:N) {
        mu[i] = foldchange(R, Nns, ep_RA[DNA_idx[i]], c[i], ep_a[IND_idx[i]], ep_i[IND_idx[i]],
                           ep_AI,n_sites);
        fc[i] ~ normal(mu[i], sigma[idx[i]])
    }
}