# ``code/stan``
This directory contains all Stan models used for parameter inference. No files
are executed in this directory. 

* ``Chure2019_DNA_binding_energy_global.stan`` \| Statistical model for
  inference of DNA binding energy for a single mutant considering *all*
  induction profiles. 
* ``Chure2019_DNA_binding_energy.stan`` \| Statistical model for inference of
  DNA binding energy for a single mutant and a single induction profile. 
* ``Chure2019_empirical_F_inference.stan`` \| Statistical model for the mean fold-change and
  standard deviation for a collection of fold-change measurements for a single
  mutant at a given repressor copy number, operator sequence, and IPTG
  concentration.
* ``Chure2019_functions.stan`` \| Contains myriad functions pertinent to this
  work. This file is called in all parameter inference models. 
* ``Chure2019_KaKi_epAI_global.stan`` \| Statistical model that infers inducer dissociation
  constants and allosteric state energy difference for a single inducer binding
  mutant considering *all* induction profiles. 
* ``Chure2019_KaKi_epAI.stan`` \| Statistical model that infers the inducer
  dissociation constants and allosteric state energy difference for an
  inducer binding mutant from a single induction profile.
* ``Chure2019_KaKi_only.stan`` \| Statistical model that infers the inducer
  dissociation constants for an inducer binding domain mutant from a single
  induction profile.