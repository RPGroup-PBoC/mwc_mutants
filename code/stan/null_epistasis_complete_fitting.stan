data {
    // Dimensional parameters
    int<lower=0> J_DNA; // Number of unique DNA mutants
    int<lower=0> J_IND; // Number of unique inducer mutants
    int<lower=0> N_DNA; // Number of DNA mutant measurements
    int<lower=0> N_IND; // Number inducer mutant measurements
    int<lower=1> N; // Total number of measurements
    int<lower=0, upper=J_DNA> DNA_idx[N_DNA]; // ID vector for DNA binding mutants
    int<lower=0, upper=J_IND> IND_idx[N_IND]; // ID vector for inducer binding mutants

    // Architectural parameters
    real ep_RA_mu; // Mean DNA binding energy for inducer mutants in kBT
    real ep_RA_sig; // Uncertainty in DNA binding energy for inducer mutants in kBT
    real R; // Number of repressors per cell for each emasurement
    real Nns; // Number of nonspecific binding sites

    // Allosteric parameters
    real ka_mu; // Mean active-repressor inducer dissociation constant in µM
    real ka_sig; // Uncertainty in active-repressor inducer dissociation constant in µM
    real ki_mu; // Mean inactive-repressor inducer dissociation constant in µM
    real ki_sig; // Uncertainty in inactive-repressor inducer dissociation constant in µM
    real ep_ai; // energetic difference between active and inactive repressor in kBT.
    int n_sites; // Nubmer of allosteric binding sites on repressor

    // Experimental parameters
    real<lower=0> DNA_c[N_DNA]; // Allosteric effector concentration for DNA binding mutant measurements
    real<lower=0> IND_c[N_IND]; // Allosteric effector concentration for DNA binding mutant measurements

    // Measurement parameters
    real DNA_foldchange[N_DNA]; // Measured fold-change for DNA binding mutants
    real IND_foldchange[N_IND]; // Measured fold-change for inducer binding mutants
}

parameters {
    // DNA binding mutant estimated parameters
    real<lower=-50, upper=0> ep_RA[J_DNA]; // DNA binding energy
    real<lower=0> DNA_sigma[J_DNA]; // Homoscdastic error in fold-change measurements

    // Inducer  binding mutant estimated parameters
    vector<lower=0, upper=1000>[J_IND] Ka; // Inducer dissoaction constant to active repressor
    vector<lower=0, upper=1000>[J_IND] Ki; // Inducer dissoaction constant to inactive repressor
    real<lower=0> IND_sigma[J_IND]; // Homoscedastic error for fold-change measurements
    
    // Define the seeded parameters.
    real DNA_ka;
    real DNA_ki;
    real IND_ep_RA;
}

transformed parameters {
    // Compute the appropriate log tranformations
    vector[J_IND] ep_a; // Active-repressor inducer binding energy
    vector[J_IND] ep_i; // Inactive-repressor inducer binding energy
    real DNA_ep_a; // Seeded value for active-repressor inducer binding energy for DNA binding mutants
    real DNA_ep_i; // Seeded value for inactive-repressor inducer binding energy for DNA binding mutants


    // Compute the log transformations
    ep_a = -log(Ka);
    ep_i = -log(Ki);
    DNA_ep_a = -log(DNA_ka);
    DNA_ep_i = -log(DNA_ki);
}

model {
    // Define the vector to compute the foldchange expected from theory.
    real mu_DNA[N_DNA]; 
    real mu_IND[N_IND];

    // Define the priors for estimated parameters. 
    DNA_sigma ~ normal(0, 1);
    IND_sigma ~ normal(0, 1);
    ep_RA ~ normal(0, 10);
    ep_a ~ normal(0, 10);
    ep_i ~ normal(0, 10);
    

    // Priors for seeded parameters.
    DNA_ka ~ normal(ka_mu, ka_sig);
    DNA_ki ~ normal(ki_mu, ki_sig);
    IND_ep_RA ~ normal(ep_RA_mu, ep_RA_sig);

    // Likelihood for DNA binding mutants.
    for (i in 1:N_DNA) {
        mu_DNA[i] = fold_change(R, Nns, ep_RA[DNA_idx[i]], DNA_c[i], DNA_ep_a,
                                DNA_ep_i, ep_ai, n_sites);
        DNA_foldchange[i] ~ normal(mu_DNA[i], DNA_sigma);
    }

    // Likelihood for inducer binding mutants.
    for (i  in 1:N_IND) {
        mu_IND[i] = fold_change(R, Nns, IND_ep_RA, IND_c[i], ep_a[IND_idx[i]],
                                ep_i[IND_idx[i]], ep_ai, n_sites);
        IND_foldchange[i] ~ normal(mu_IND[i], IND_sigma);
    }

}
