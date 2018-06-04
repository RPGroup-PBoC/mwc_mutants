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

  // Dimensional parameters
  int N; // Total number of data points.

  // Architectural parameters
  vector<lower=0>[N] R; // Number of repressors
  real<lower=0> Nns; // Number of nonspecific binding sites

  // Allosteric parameters 
  vector<lower=0>[N] c; // Effector concentration.
  real epR; // Binding energy in kBT. 
  real ep_ai; // Allosteric energy difference 
  int<lower=1> n_sites; // Number of allosteric sites.  

  // Observed parameters.
  vector<lower=-0.1, upper=1.2>[N] fc;
  }

parameters {
  real<lower=0> ka; // Active repressor inducer dissociation constant
  real<lower=0> ki; // Inactive repressor inducer dissociation constant
  real<lower=0> sigma; //  Homoscedastic error
}

transformed parameters {
  real<lower=-7, upper=7>ep_a; 
  real<lower=-7, upper=7>ep_i;
  ep_a = log(ka);
  ep_i = log(ki);
}

model {
  vector[N] mu;

  // Define the priors. 
  sigma ~ normal(0, 1);
  ep_a ~ normal(0, 10);
  ep_i ~ normal(0, 10);

  for (i in 1:N) {
    mu[i] = fold_change(R[i], Nns, epR, c[i], ep_a, ep_i, ep_ai, n_sites);
  }

  // Define the likelihood
  fc ~ normal(mu, sigma);
}