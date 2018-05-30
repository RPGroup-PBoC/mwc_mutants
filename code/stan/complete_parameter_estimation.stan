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
  // Hierarchical parameters which provide dimensionality information
  int<lower=0> J_DNA; // Number of distinct DNA binding mutants
  int<lower=0> J_IND; // Number of distinct inducer binding mutants
  int<lower=0> J_DBL; // Number of distinct double mutants
  int<lower=0> N_DNA; // Total number of data points for DNA measurements
  int<lower=0> N_IND; // Total number of data points for IND measurements
  int<lower=0> N_DBL; // Total number of data points for DBL measurements

  // Vectors of trial identifiers
  int<lower=0, upper=J_DNA> idx_DNA[N_DNA];
  int<lower=0, upper=J_IND> idx_IND[N_IND];
  int<lower=0, upper=J_DBL> idx_DBL[N_DBL];
 
  // Generate vectors for repressor counts
  vector<lower=0>[N_DNA] R_DNA; 
  vector<lower=0>[N_IND] R_IND; 
  vector<lower=0>[N_DBL] R_DBL; 

  // Define wild-type constants
  real<lower=0> wt_ka; // In units of µM
  real<lower=0> wt_ki; // In units of µM
  real wt_epR; // In units of kBT

  // Define thermodynamic parameters
  real ep_ai;
  real<lower=0> Nns; 
  int<lower=0> n_sites;
  vector<lower=0>[N_DNA] c_DNA; // IPTG concentration for DNA mutants. 
  vector<lower=0>[N_IND] c_IND; // IPTG concentration for DNA mutants. 
  vector<lower=0>[N_DBL] c_DBL; // IPTG concentration for DNA mutants. 

  // Define measurement parameters
  real fc_DNA[N_DNA]; // Fold-change measurements for DNA mutants
  real fc_IND[N_IND]; // Fold-change measurements for IND mutants
  real fc_DBL[N_DBL]; // Fold-change measurements for DBL mutants  
  }

parameters { 
  // Define fitting parameters for DNA and inducer mutants
  real DNA_epR[J_DNA]; 
  real<lower=0> IND_ka[J_IND];
  real<lower=0> IND_ki[J_IND];

  // Define fitting parameters for double mutants
  real DBL_epR[J_DBL];
  real<lower=0> DBL_ka[J_DBL];
  real<lower=0> DBL_ki[J_DBL];

  // Homoscedastic error definition
  vector<lower=0>[J_DNA] DNA_sigma;
  vector<lower=0>[J_IND] IND_sigma;
  vector<lower=0>[J_DBL] DBL_sigma;

  }

transformed parameters{
  // Define log_transformations
 real log_ka_DNA; 
 real log_ki_DNA;
 real log_ka_IND[J_IND];
 real log_ki_IND[J_IND];
 real log_ka_DBL[J_DBL];
 real log_ki_DBL[J_DBL];

  // Compute tranformations
  log_ka_DNA = log(wt_ka);
  log_ki_DNA = log(wt_ki);
  log_ka_IND = log(IND_ka);
  log_ki_IND = log(IND_ki);
  log_ka_DBL = log(DBL_ka);
  log_ki_DBL = log(DBL_ki);
  }

model {
  // Assign vectors for calculation of fold-change.
  vector[N_DNA] DNA_mu;
  vector[N_IND] IND_mu;
  vector[N_DBL] DBL_mu;

  // Thermodynamic model parameter prior definition
  DNA_epR ~ normal(0, 10);
  DBL_epR ~ normal(0, 10);
  IND_ka ~ normal(0, 100);
  IND_ki ~ normal(0, 100);
  DBL_ka ~ normal(0, 100);
  DBL_ki ~ normal(0, 100);

  // Homoscedastic error prior definition
  DNA_sigma ~ normal(0, 1);
  IND_sigma ~ normal(0, 1);
  DBL_sigma ~ normal(0, 1);

  // DNA binding mutant likelihood
  for (i in 1:N_DNA){
    DNA_mu[i] = fold_change(R_DNA[i], Nns, DNA_epR[idx_DNA[i]], c_DNA[i], log_ka_DNA,
                        log_ki_DNA, ep_ai, n_sites);
    fc_DNA[i] ~ normal(DNA_mu[i], DNA_sigma[idx_DNA[i]]);
    }

  // Inducer binding likelihood
  for (i in 1:N_IND) {
    IND_mu[i] = fold_change(R_IND[i], Nns, wt_epR, c_IND[i], -log_ka_IND[idx_IND[i]], 
                            -log_ki_IND[idx_IND[i]], ep_ai, n_sites);
    fc_IND[i] ~ normal(IND_mu[i], IND_sigma[idx_IND[i]]);
    }

  // DBL binding likelihood
  for (i in 1:N_DBL) {
    DBL_mu[i] = fold_change(R_DBL[i], Nns, DBL_epR[idx_DBL[i]], c_DBL[i], -log_ka_DBL[idx_DBL[i]], 
                        -log_ki_DBL[idx_DBL[i]], ep_ai, n_sites);
    fc_DBL ~ normal(DBL_mu[i], DBL_sigma[idx_DBL[i]]);
    }
  }