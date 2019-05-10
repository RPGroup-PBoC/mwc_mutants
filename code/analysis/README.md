# ``code/analysis``

This repository contains all code that performs analysis routines. The files are
intended to be run from this directory and have relative paths to data and Stan models.

* `Chure2019_DNA_prior_predictive_samples.py` \| Generates samples for prior
  predictive checks of the DNA binding energy inference for the DNA binding
  domain mutants. 
* `Chure2019_DNA_sbc_samples.py` \| Generates statistics of simulation based calibration
  and samples for DNA binding energy inference for the DNA binding
  domain mutants. 
* `Chure2019_empirical_F_inference_wt.py` \| Infers the mean free energy from a
  collection of fold-change measurements of wild-type data from Razo-Mejia et
  al. 2018.
* `Chure2019_empirical_F_limits.py` \| Generates simulated data and MCMC samples
  for exploring the inferential limits of the free energy from fold-change data.
* `Chure2019_empirical_F_prior_predictive_checks.py` \| Generate prior
  predictive checks for inference of the free energy from fold-change
  measurements.
* `Chure2019_empirical_F_sbc_samples.py` \| Generates samples and statistics for
  simulation based calibration of empirical free energy inference. 
* `Chure2019_IND_posterior_predictive_checks.py` \| Executes posterior
  predictive checks for both inferential models for the inducer binding domain
  mutants. 
* `Chure2019_IND_prior_predictive_samples.py` \| Generates prior predictive
  samples for both hypotheses of inducer binding domain mutants. 
* `Chure2019_IND_sbc_samples.py` \| Generates statistics and samples for
  simulation based calibration of both inferential models for inducer binding
  domain mutants. 
* `Chure2019_parameter_inference.py` \| Performs all parameter inference *on
  single induction profiles* for each mutant. This script calls the `stan`
  models which live in another subdirectory. 
* `Chure2019_pooled_parameter_inference.py` \| Performs parameter inference for
  each mutant considering *all* induction profiles. This estimates the DNA
  binding energy for DNA binding mutants and allosteric parameters for inducer
  binding domain mutants.