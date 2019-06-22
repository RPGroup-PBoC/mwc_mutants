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
    real repression(real pact, real R, real Nns, real ep_r, real offset) {
        return 1 + pact * (R / Nns) * exp(-ep_r) + (1 - pact) * (R / Nns) *
        exp(-(ep_r + offset));
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
                    real ep_ai, int n_sites, real offset) {
        // Compute the various componenets piecewise for simplicity.
        real pact;
        real rep;
        pact = prob_act(c, ep_a, ep_i, ep_ai, n_sites);
        rep = repression(pact, R, Nns, ep_r, offset);
        return rep^-1;
        }
      }
data {
    // Dimensional information
    int<lower=0> N; // total number of data points
    int<lower=0> J; // Number of unique operators
    int<lower=0, upper=J> idx[N]; // Identification vector for operator

    // Preset constants
    real ep_ai; // Energy difference between states
    real<lower=0> n_ns; // Number of nonspecific binding sites
    int<lower=1> n; // Number of allosteric sites.

   // Observed data
   vector[N] fc; // Observed fold-change
   real<lower=0> c[N]; // Inducer concentration
   real<lower=0> R[N]; // Repressors per cell
}

parameters {
    // Define the lac specific parameters
    real ep_a; // Log transform of Ka
    real ep_i; // Log transform of Ki
    real sigma; // Homoscedastic error
    vector[J] ep_r; // DNA binding energy
}

transformed parameters {
    real ka = exp(ep_a);
    real ki = exp(ep_i);
}

model {
    // Compute the theoretical mean foldchange
    vector[N] mu;
    for (i in 1:N) { 
        mu[i] = fold_change(R[i], n_ns, ep_r[idx[i]], c[i], ep_a, ep_i, ep_ai,
        n, log(1000));
    }

    // Set the priors
    ep_a ~ normal(0, 3);
    ep_i ~ normal(0, 3);
    ep_r ~ normal(-12, 6);
    sigma ~ normal(0, 0.1);

    // Evaluate the likelihood
    fc ~ normal(mu, sigma);
}