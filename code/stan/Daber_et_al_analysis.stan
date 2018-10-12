/* 
* Analysis model for Daber et al. 
* ---------------------------------------------------------------------
* Author: Griffin Chure
* License: MIT
* 
* Description
* ---------------------------------------------------------------------
* This model samples the posterior probability distribution for the
* parameters Ka, Ki, and R/Kdna.
*/

data {
  // Dimensional parameters
  int<lower=1> J_DNA; // Number of unique DNA binding domain mutants
  int<lower=1> J_IND; // Number of unique inducer binding domain mutants
  int<lower=1> N_DNA; // Number of measurements of DNA binding mutants
  int<lower=1> N_IND; // Number of measurements of inducer binding mutants
  int<lower=1> N_WT; // Number of measurements of wild-type strain
  
  // Identification vectors
  int<lower=1, upper=J_DNA> DNA_idx[N_DNA];
  int<lower=1, upper=J_IND> IND_idx[N_IND];
      
  // Known parameters
  int n_sites; // Number of allosteric binding sites of the repressors
  real ep_AI; // Energetic difference between active adn inactive allosteric states in k_BT
  
  // IPTG concentrations for each sample in µM
  real<lower=0> c_DNA[N_DNA]; 
  real<lower=0> c_IND[N_IND];
  real<lower=0> c_WT[N_WT];
  
  // Observed fold-change.
  real DNA_fc[N_DNA];
  real IND_fc[N_IND];
  real WT_fc[N_WT];
}

parameters {
    // Wild-type parameters.
    real<lower=0> R_Kdna_wt;
    real<lower=0> Ka_wt;
    real<lower=0> Ki_wt;
    
    // DNA binding parameters
    real<lower=0> R_Kdna_mut[J_DNA];
    
    // Inducer binding parameters
    real<lower=0> Ka_mut[J_IND];
    real<lower=0> Ki_mut[J_IND];
    
    // Homoscedastic error
    real<lower=0> sigma_wt;
    real<lower=0> sigma_dna[J_DNA];
    real<lower=0> sigma_ind[J_IND];
}

transformed parameters {
    // Perform a log transformation of the inducer binding constants for better sampling. 
    real log_ka_wt;
    real log_ki_wt;
    real log_ka_mut[J_IND];
    real log_ki_mut[J_IND];
    log_ka_wt = log(Ka_wt);
    log_ki_wt = log(Ki_wt);
    for (i in 1:J_IND) {
        log_ka_mut[i] = log(Ka_mut[i]);
        log_ki_mut[i] = log(Ki_mut[i]);
    }
}

model {
    // Define a vector for calculation of fold_change
    vector[N_WT] mu_wt;
    vector[N_DNA] mu_dna;
    vector[N_IND] mu_ind;
    
    // Define placeholders to compute the p_act
    real p_act_wt;
    real p_act_dna;
    real p_act_ind;
    
    // Define the various priors. 
    R_Kdna_wt ~ lognormal(0, 3);
    R_Kdna_mut ~ lognormal(0, 3);
    log_ka_wt ~ normal(0, 1);
    log_ki_wt ~ normal(0, 1);
    log_ka_mut ~ normal(0, 1);
    log_ki_mut ~ normal(0, 1);
    sigma_wt ~ normal(0, 1);
    sigma_dna ~ normal(0, 1);
    sigma_ind ~ normal(0, 1);
    
    // Define the likelihood for the wild-type samples. 
    for (i in 1:N_WT) {
       // Compute the probability of an active repressor at the concentration. 
       p_act_wt = (1 + c_WT[i] * exp(-log_ka_wt))^n_sites / ((1 + c_WT[i] * exp(-log_ka_wt))^n_sites + exp(-ep_AI) * (1 + c_WT[i] * exp(-log_ki_wt))^n_sites);
       mu_wt[i] = (1 + p_act_wt * R_Kdna_wt)^-1;
       WT_fc[i] ~ normal(mu_wt[i], sigma_wt);
    }
    
    // Define the likelihood for the DNA binding mutants. 
    for (i in 1:N_DNA) {
        p_act_dna = (1 + c_DNA[i] * exp(-log_ka_wt))^n_sites / ((1 + c_DNA[i] * exp(-log_ka_wt))^n_sites + exp(-ep_AI) * (1 + c_DNA[i] * exp(-log_ki_wt))^n_sites);
        mu_dna[i] = (1 + p_act_dna * R_Kdna_mut[DNA_idx[i]])^-1;
        DNA_fc[i] ~ normal(mu_dna[i], sigma_dna[DNA_idx[i]]);
    }
   
   // Define the likelihood for the inducer binding mutants.
   for (i in 1:N_IND) {
       p_act_ind = (1 + c_IND[i] * exp(-log_ka_mut[IND_idx[i]]))^n_sites / ((1 + c_IND[i] * exp(-log_ka_mut[IND_idx[i]]))^n_sites + exp(-ep_AI) * (1 + c_IND[i] * exp(-log_ki_mut[IND_idx[i]]))^n_sites);
       mu_ind[i] = (1 + p_act_ind * R_Kdna_wt)^-1;
       IND_fc[i] ~ normal(mu_ind[i] , sigma_ind[IND_idx[i]]);
   }
}