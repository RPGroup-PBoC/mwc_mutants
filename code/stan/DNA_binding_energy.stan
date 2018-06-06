functions {


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
    real fold_change(real R, real Nns, real ep_r, real ep_ai) {
        // Compute the various componenets piecewise for simplicity.
        real pact;
        real rep;
        pact = (1 + exp(-ep_ai))^-1;
        rep = repression(pact, R, Nns, ep_r);
        return rep^-1;
        }
      }

data { 
  // Dimensional parameters
  int J; // Number of unique mutants
  int N; // Total number of data points.
  int idx[N]; // Trial idx

  // Architectural parameters
  vector<lower=0>[N] R; // Number of repressors
  real<lower=0> Nns; // Number of nonspecific binding sites

  // Allosteric parameters 
  real ep_ai; // Allosteric energy difference 

  // Observed parameters.
  vector<lower=-0.1, upper=1.2>[N] fc;
  }

parameters {
  real epR[J]; // DNA binding energy in units of kBT
  real<lower=0> sigma[J]; //  Homoscedastic error
}

transformed parameters {
  vector[N] log_fc;
  log_fc = log(fc);
}

model {
  vector[N] mu;

  // Define the priors. 
  sigma ~ normal(0, 1);
  epR ~ normal(0, 10);

  for (i in 1:N) {
    mu[i] = fold_change(R[i], Nns, epR[idx[i]], ep_ai);
    log_fc ~ normal(log(mu[i]), log(sigma[idx[i]]));
  }
}