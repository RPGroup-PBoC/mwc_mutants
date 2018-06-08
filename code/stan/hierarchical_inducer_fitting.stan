functions {
    /**
    * Compute the probability of a repressor being active given an inducer
    * concentration c.
    *
    * @param c Concentration of allosteric effector.
    * @param ep_a Log transform of effector dissociation constant from active
    *        repressor, Ka, in kBT.
    * @param ep_i Log transform of effector dissociation constant from inactive
    *        repressor, Ki, in kBT.
    * @param ep_ai Energy difference between the active and inactive state of
    *        the repressor in kBT.
    * @param n_sites The number of allosterically dependent sites.
    * @return prob_act The probability of a repressor being active with the
    *         given parameters.
    **/
    real prob_act(real c, real ep_a, real ep_i, real ep_ai, int n_sites) {
        // Calculate the relevant components piecewise for simplicity.
        real numerator;
        real denominator;
        numerator = (1 + c * exp(-ep_a))^n_sites;
        denominator = numerator + exp(-ep_ai) * (1 + c * exp(-ep_i))^n_sites;
        return numerator / denominator;}

    /**
    * Compute the level of repression in a simple repression architecture.
    *
    * @param pact The probability of an active repressor.
    * @param R The number of repressors per cell.
    * @param Nns The number of nonspecific binding sites.
    * @param ep_r The binding energy of the repressor to the DNA in kBT.
    * @return repression The level of repression given these parameters.
    **/
    real repression(real pact, real R, real Nns, real ep_r) {
        return 1 + pact * (R / Nns) * exp(-ep_r);
      }

    /**
    * Calculate the fold-change in gene expression.
    *
    * @param R The number of repressors per cell
    * @param Nns The number of nonspecific repressor binding sites.
    * @param ep_r The  binding energy of the repressor to the DNA in kBT.
    * @param c The concentration of allosteric effector.
    * @param ep_a The log transform of the effector dissociation constant from
    *        the active repressor, Ka, in kBT.
    * @param ep_i The log tranform of the effector dissociation constant from
    *        the active repressor, Ki, in kBT.
    * @param ep_ai The energetic difference between the active and inactive
    *        states of the repressor in kBT.
    * @param n_sites The number of allostericaly dependent effector binding
    *        sites.
    **/
    real fold_change(real R, real Nns, real ep_r, real c, real ep_a, real ep_i,
                    real ep_ai, int n_sites) {
        // Compute the various componenets piecewise for simplicity.
        real pact;
        real rep;
        pact = prob_act(c, ep_a, ep_i, ep_ai, n_sites);
        rep = repression(pact, R, Nns, ep_r);
        return rep^-1;
        }
      }

data {
    // Dimensional parameters.
    int<lower=1> J; // Number of unique mutants
    int<lower=0> N_leak; // Number of leakiness points
    int<lower=0> N; // Total number of points
    int<lower=0, upper=J> idx_leak[N_leak]; // Identification vector for leakiness
    int<lower=0, upper=J> idx[N]; // Identification vector for titration curve

    // Architectural parameters
    real<lower=0> R; // Number of repressors
    real<lower=0> Nns; // Number of nonspecific binding sites
    real epR; 

    // Allosteric parameters
    int n_sites; // Number of allosteric binding sites
    vector<lower=0>[N] c; 
    // Measured parameters. 
    vector<lower=-0.1, upper=1.2>[N_leak] fc_leak; // Leakiness fold-change measurements.
    vector<lower=-0.1, upper=1.2>[N] fc; // All fold-change measurements
}

parameters {
    vector[J] ep_ai;
    vector<lower=0>[J] ka;
    vector<lower=0>[J] ki;
    vector<lower=0>[J] sigma;
}

transformed parameters {
    // Peform log transformation of ka, ki, and fc for leakiness
    vector[J] ep_a;
    vector[J] ep_i;
    ep_a = log(ka);
    ep_i = log(ki);
}

model {
    // Set up vectorws
    vector[N_leak] mu_leak;
    vector[N] mu;

    // Set the priors
    ep_ai ~ normal(0, 1);
    ep_a ~ normal(0, 1);
    ep_i ~ normal(0, 1);
    sigma ~ normal(0, 1);

    // Compute the likelihood for ep_ai.
    for (i in 1:N_leak) {
        mu_leak[i] = fold_change(R, Nns, epR, 0, ep_a[idx_leak[i]], ep_i[idx_leak[i]], 
                                ep_ai[idx_leak[i]], n_sites);
        fc_leak[i] ~ normal(mu_leak[i], sigma);
    }

    // Compute the likelihood for ka, ki
    for (i in 1:N) {
        mu[i] = fold_change(R, Nns, epR, c[i], ep_a[idx[i]], ep_i[idx[i]], ep_ai[idx[i]], n_sites);
        fc[i] ~ normal(mu[i], sigma[idx[i]]);
    }

}