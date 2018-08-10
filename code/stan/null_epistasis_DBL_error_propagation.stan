/*
* Double mutant error propgration assuming no epistasis
* -----------------------------------------------------
* Author: Griffin Chure
* Date: 08 August 2018
* License: CC-BY 4.0
*
* Description
* -----------
* This model estimates the DNA binding energy as well as the inducer dissociation constants
* for the double mutants given measurements of the fold-change. The parameter estimates are 
* performed using informative priors based on the mode and HPD for parameter estimates from 
* the single mutants. This model assumes no epistasis between mutants.
*/
data {
    // Dimensional parameters
    int<lower=1> J_DNA; // Number of unique DNA-binding domain mutants
    int<lower=1> J_IND; // Number of unique inducer-binding domain mutants
    int<lower=1> N; // Total number of data points.
    int<lower=1, upper=J_DNA> DNA_idx[N]; // ID vector for DNA binding domain mutants
    int<lower=1, upper=J_IND> IND_idx[N]; // ID vector for inducer-binding domain mutants

    // Architectural parameters
    real ep_RA_mode[J_DNA]; // Best-fit values for DNA mutants
    real ep_RA_lower[J_DNA]; // Upper bound for DNA binding energy
    real ep_RA_upper[J_DNA]; // Lower bound for DNA binding energy
    real R[N]; // Repressor copy number
    real Nns; // Number of nonspecific binding sites

    // Allosteric parameters (dissociation constants in units of µM)
    real ka_mode[J_IND]; // Best-fit values for inducer dissocataion constant to active repressor
    real ka_lower[J_IND]; // Lower bound for for inducer dissociation constant to active repressor 
    real ka_upper[J_IND]; // Upper bound for for inducer dissociation constant to active repressor 
    real ki_mode[J_IND]; // Best-fit values for inducer dissocataion constant to inactive repressor
    real ki_lower[J_IND]; // Lower bound for for inducer dissociation constant to inactive repressor 
    real ki_upper[J_IND]; // Upper bound for for inducer dissociation constant to inactive repressor 
    int n_sites; // Number of allosteric effector binding sites 
    real c[N]; // Allosteric effector concentration in µM

    // Measured parameters
    vector[N] foldchange; // Measured fold-change in gene expression
}

parameters {

}