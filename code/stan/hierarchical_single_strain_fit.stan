data {
    // Information about data dimensions.
    int<lower=1> J; // Number of unique datasets
    int<lower=1> N; // Total number 
    int<lower=1, upper=J> trial[N];

    // Architectural parameters
    vector<lower=0>[N] R; // Number of repressors.
    real<lower=1> n_ns; // Number of nonspecific binding sites

    // Allosteric parameters
    real ep_ai; // Allosteric energy difference. 
    int<lower=1> n_sites; // Number of allosteric sites.    
    vector<lower=0>[N] c; // Effector concentration.
    
    // Measureables
    vector[N] fc;
}

parameters {
    real ep_R[J]; // DNA binding energy
    real<lower=1E-12> ka[J]; // Active repressor inducer dissociation constant. 
    real<lower=1E-12> ki[J]; // Inactive repressor inducer dissociation constant.
    real<lower=1E-9> sigma[J]; // Homoscedastic error
}

transformed parameters {
    real ka_tilde[J];
    real ki_tilde[J];
    ka_tilde = log(ka);
    ki_tilde = log(ki);
}

model {
    vector[N] mu;
    ep_R ~ normal(0, 100);
    ka ~ normal(0, 100);
    ki ~ normal(0, 100);
    sigma ~ normal(0, 10);

    for (i in 1:N) {
        mu[i] = fold_change(R[i], n_ns, ep_R[trial[i]], c[i], ka_tilde[trial[i]],
                            ki_tilde[trial[i]], ep_ai, n_sites);
        fc[i] ~ normal(mu[i], sigma[trial[i]]);
    }
}