data {
    // Heirarchical data parameters
    int<lower=0> J; // Number of trials (mutants).
    int<lower=0> N; // total number of measurements.
    int<lower=1, upper=J> trial[N]; // trial ID for each measurement

    // Simple repression parameters.
    real<lower=0> R[N]; // Number of repressors
    real<lower=0> n_ns; // Number of nonspecific binding sites 

    // Allosteric parameters
    real<lower=0> c[N]; // effector concentration in units of µM
    real ka; // in units of µM
    real ki; // in units of µM
    real ep_AI; // Allosteric energy difference in kBT
    int<lower=1> n_sites; // Number of allosteric sites

    // Measurement parameters
    real<lower=-0.2, upper=1.3> fc[N]; // Experimentall observed fold-change
 }

parameters {
    real ep_R[J]; // Repressor-DNA binding energy in kBT
    real<lower=1E-9> sigma[J]; // Homoscedastic error. Assumed constant for all.
}

model {
    // Compute the expected value of the fold-change
    vector[N] mu;

    // Set the priors to be weakly informative
    ep_R ~ normal(0, 10);
    sigma ~ normal(0, 10);

    for (i in 1:N) {
        mu[i] = fold_change(R[i], n_ns, ep_R[trial[i]], c[i], 
                            ka, ki, ep_AI, n_sites);
        // Define the likelihood.
         fc[i] ~ normal(mu[i], sigma[trial[i]]);

    }

 }

