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
    // Define dimensionality
    int J; // Number of unique double mutants.
    int N; // Total number of measurements
    real<lower=0, upper=J> idx[N]; // Vector of mutant identifiers. 

    // Define the constant parameters
    real<lower=0> R; 
    real<lower=0> Nns;
    real ep_ai;
    real n_sites;
    real<lower=0> c[N]

    // Define the parameters to establish priors. 
    real ka_mu[J];
    real ki_mu[J];
    real epR_mu[J];
    real ka_sig[J];
    real ki_sig[J];
    real epR_sig[J];

    // Define the measured parameter
    real fc[N]; 
}

parameters {
    # Define the fitting parameters
    real<lower=0> ka[J];
    real<lower=0> ki[J];
    real epR[J];
    real<lower=0> sigma[J]; 
}

model {
    // Define a vector in which to compute the theoretical fold-change.
    vector[N] mu;

    // Define the priors using the informative values.
    ka ~ normal(ka_mu, ka_sig);
    ki ~ normal(ki_mu, ki_sig);
    epR ~ normal(epR_mu, epR_sig);

    // Define homoscedastic error prior
    sigma ~ normal(0, 10);

    // Compute the fold-change.
    for (i in 1:N) {
        mu[i] = fold_change(R, Nns, epR[idx[i]], c[i], log(ka[idx[i]]), log(ki[idx[i]]),
                            ep_ai, n_sites);
        fc ~ normal(mu[i], sigma[idx[i]]);
    }
}
