data {
    // Dimension details.
    int<lower=1> J; // Number of unique sets
    int<lower=1> N; // Number of measurements
    int<lower=1, upper=J> trial[N]; // Vector with trial identifier

    // Allosteric details
    real ep_AI; // Allosteric energy difference
    real ka; // Active repressor inducer dissociation constant (µM)
    real<lower=1E-9> ka_sigma; // Std. dev from Gaussian approximation of ka
    real<lower=0> ki; // Inactive repressor inducer dissociation constant (µM)
    real ki_sigma; // Std. dev from Gaussian approximation of ki
    real c[N]; // Concentration of effector
    int<lower=1> n_sites; // Number of allosteric sites

    // Architectural details
    real<lower=0> R[N]; // Number of repressors per cell
    real<lower=1> n_ns; // Number of nonspecific binding sites
    real fc[N]; // measured output
    }

parameters {
    real ep_R[J];  // Repressor binding energy
    real<lower=1E-9> sigma[J]; // Homoscedastic error
    }

model {
    // Compute the expected foldchange and likelihood.
    vector[N] mu; 

    real ka_dist;
    real ki_dist;
    ka_dist ~ normal(ka, ka_sigma);
    ki_dist ~ normal(ki, ki_sigma);
    
    // Set the priors
    ep_R ~  normal(0, 100); // Weakly informative
    sigma ~  normal(0, 100);// Weakly informative

    for (i in 1:N) {
        mu[i] = fold_change(R[i], n_ns, ep_R[trial[i]], c[i], -log(ka_dist), -log(ki_dist),
                            ep_AI, n_sites);
        fc[i] ~ normal(mu[i], sigma[trial[i]]);
        }
    }