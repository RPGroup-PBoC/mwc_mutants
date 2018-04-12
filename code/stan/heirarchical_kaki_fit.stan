data {
    // Dimensional parameters
    int<lower=1> J; // Number of unique mutants
    int<lower=1> N; // Number of data points.
    int<lower=1, upper=J> trial[N]; // Vector of trial identifiers.

    // Architectural parameters.
    real R; // Number of repressors 
    real epR; // DNA binding energy
    real<lower=1> n_ns; // Number of nonspecific binding sites

    // Allosteric parameters
    real c[N]; // Effector concentration
    real ep_ai; // Allosteric energy difference
    int<lower=1> n_sites; // Number of allosteric sites.

    // Measured parameters
    real fc[N]; // observed fold-change.
}

parameters {
    real ka[J]; // Inducer dissociation constant to active repressor
    real ki[J]; // Inducer dissociation constant to inactive repressor
    real<lower=0> sigma[J];  // Homoscedastic error

}

model {
    vector[N] mu; 

    // Priors
    ka ~ normal(0, 100);
    ki ~ normal(0, 100);
    sigma ~ normal(0, 100);

    // Likelihood
    for (i in 1:N) {
        mu[i] = fold_change(R, n_ns, epR, c[i], ka[trial[i]], ki[trial[i]], ep_ai, n_sites);
        fc[i] ~ normal(mu[i], sigma[trial[i]]);
    }
}