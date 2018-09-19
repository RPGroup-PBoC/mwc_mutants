# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import pystan
sys.path.insert(0, '../../')
import mut.bayes
import mut.stats
import tqdm
constants = mut.thermo.load_constants()

# Define the figure constants
FIG_NO = 5
OPERATOR = 'O2'
RBS = 'RBS1027'
n_draws = int(1E5)

# Load the raw data
data = pd.read_csv('../../data/csv/compiled_data.csv')
muts = data[data['class']=='DBL']['mutant'].unique()

# Load the samples
epRA_samples = pd.read_csv('../../data/csv/Fig2_O2_DNA_binding_energy_samples.csv')
kaki_samples = pd.read_csv('../../data/csv/Fig3_O2_KaKi_only_samples.csv')
kaki_epAI_samples = pd.read_csv('../../data/csv/Fig3_O2_KaKi_epAI_samples.csv')

# Identify the double mutants. 
dfs = []
DNA_muts = []
IND_muts = []
for m in muts:
    dna, ind = m.split('-')
    DNA_muts.append(dna)
    IND_muts.append(ind)
DNA_muts = list(set(DNA_muts))
IND_muts = list(set(IND_muts))

# Iterate through and compute the saturation fold-change. 
def delta_bohr(mut_constants, operator, wt_constants):
    DNA_domain = wt_constants[operator] - mut_constants['epRA']
    IND_numer = 1 + np.exp(-mut_constants['ep_AI']) *\
                   (1 + (mut_constants['Ka']/mut_constants['Ki']))**2
    DNA_denom = 1 + np.exp(-wt_constants['ep_AI']) *\
                    (1 + (wt_constants['Ka']/wt_constants['Ki']))**2
    return DNA_domain - np.log(DNA_denom) + np.log(IND_numer)

bohr_df = pd.DataFrame([])
for i, dna in enumerate(DNA_muts):
    epRA_samps = epRA_samples['ep_RA.{}.{}'.format(dna, constants[RBS])].sample(n_draws, replace=True).values
    for j, ind in enumerate(IND_muts):
        if ind == 'Q294K':
            _kaki_samps = kaki_epAI_samples.sample(n_draws, replace=True)
            epAI_samps = _kaki_samps['ep_AI.{}'.format(ind)].values
        else:
            _kaki_samps = kaki_epAI_samples.sample(n_draws, replace=True)
            epAI_samps = np.ones(n_draws) * constants['ep_AI']
        ka_draws = _kaki_samps['Ka.{}'.format(ind)].values
        ki_draws = _kaki_samps['Ki.{}'.format(ind)].values
        
        # Compute the delta bohr for each draw
        _muts = dict(epRA=epRA_samps, ep_AI=epAI_samps, Ka=ka_draws, Ki=ki_draws)
        dBohr = delta_bohr(_muts, 'O2', constants)
        cred_region = mut.stats.compute_hpd(dBohr, 0.95)
        bohr_df = bohr_df.append({'hpd_min':cred_region[0], 'hpd_max':cred_region[1], 'mutant':'{}-{}'.format(dna, ind)},
                                ignore_index=True)
bohr_df.to_csv('../../data/csv/Fig{}_{}_DBL_predicted_deltaBohr.csv'.format(FIG_NO, OPERATOR))
        
        
# Compute the measured (through fitting) deltaBohr        

        
        
        
        