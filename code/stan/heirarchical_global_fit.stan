data {
    // Dimensional information
    int<lower=1> J; // Number of independent elements
    int<lower=1> N; // Number of distinct measurements
    int<lower=1, upper=J> trial[N]; // Vector of trial ids

    // Architectural information
    vector[N] R; // Repressors per cell
    real<lower=1> n_ns; // Number of nonspecific binding sites
    
    // Allosteric information
    real ep_ai; // Allosteric energy difference
    int<lower=1> n_sites; // Number of allosteric binding sites
    vector[N] c; // IPTG concentration
    vector[N] fc; // Measured fold-change
}

parameters {
    vector[J] ep_R; // DNA binding energy in k_BT
    vector<lower=0>[J] ka; // Active repressor inducer dissociation constant
    vector<lower=0>[J] ki; // Inactive repressor inducer dissociation constant
    vector<lower=0>[J] sigma; // Homoscedastic error
}

transformed parameters {
    vector[J] ka_tilde;
    vector[J] ki_tilde;
    ka_tilde = log(ka);
    ki_tilde = log(ki);
}

model {
    vector[N] mu;
    ep_R ~ normal(0, 100);
    ka ~ normal(0, 100);
    ki ~ normal(0, 100);
    sigma ~ normal(0, 100);

    for (i in 1:N) {
        mu[i] = fold_change(R[i], n_ns, ep_R[trial[i]], c[i], 
                         ka_tilde[trial[i]], ki_tilde[trial[i]], ep_ai, n_sites);
        fc[i] ~ normal(mu[i], sigma[trial[i]]);
    }
}