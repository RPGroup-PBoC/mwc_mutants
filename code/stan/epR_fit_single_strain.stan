data {
    // Experimental parameters
    int<lower=0> N; // Number of measurements

    // Architectural parameters
    real<lower=0> R; // Number of repressors
    real<lower=1> n_ns; // Number of nonspecific binding sites.

    // Allosteric parameters.
    vector<lower=0>[N] c; // Effector concentration
    real<lower=0> ka; // Active repressor inducer dissociation constant
    real<lower=0> ki; // Inactive repressor inducer dissociation constant
    real ep_ai; // Allosteric energy difference
    int<lower=1> n_sites; // Number of allosteric binding sites.

    // Observed paraemeters
    vector[N] fc; // Observed fold-change.
}

parameters {
    real ep_R; // DNA binding energy in k_BT. 
    real<lower=1E-9> sigma; // Homoscedastic error
}

transformed parameters {
    // Transformation of ka/ki into log space
    real ep_a;
    real ep_i;
    ep_a = log(ka);
    ep_i = log(ki);
}

model {
    // Compute the expected value
    vector[N] mu;
    mu = fold_change(R, n_ns, ep_R, c, ep_a, ep_i, ep_ai, n_sites);

    // Define the priors
    ep_R ~ normal(0, 100);
    sigma ~ normal(0 10);

    // Define the likelihood
    fc ~ normal(mu, sigma)
}