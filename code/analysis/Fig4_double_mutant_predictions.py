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

# Constants for file naming. 
FIG_NO = 4
OPERATOR = 'O2'
RBS = 'RBS1027'

# %% Load the sampling traces. 
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[data['class']=='DBL']
epRA_samples = pd.read_csv('../../data/csv/Fig2_{}_DNA_binding_energy_samples.csv'.format(OPERATOR))
epRA_stats = pd.read_csv('../../data/csv/Fig2_{}_DNA_binding_energy_stats.csv'.format(OPERATOR))
kaki_samples = pd.read_csv('../../data/csv/Fig3_{}_KaKi_only_samples.csv'.format(OPERATOR))
kaki_stats = pd.read_csv('../../data/csv/Fig3_{}_KaKi_only_stats.csv'.format(OPERATOR))
kaki_epAI_samples = pd.read_csv('../../data/csv/Fig3_{}_KaKi_epAI_samples.csv'.format(OPERATOR))
kaki_epAI_stats = pd.read_csv('../../data/csv/Fig3_{}_KaKi_epAI_stats.csv'.format(OPERATOR))

# Define the sampling constants. 
n_draws = int(1E5)

# Get the mutants
muts = data['mutant'].unique()
IND_muts = []
DNA_muts = []
for m in muts:
    dna, ind = m.split('-')
    IND_muts.append(ind)
    DNA_muts.append(dna)
DNA_muts = list(set(DNA_muts))
IND_muts = list(set(IND_muts))

# Loop through mutants and isolate parameters.
dfs = []
c_range = np.logspace(-2, 4, 200)
for i, dna in enumerate(tqdm.tqdm(DNA_muts)):
    epRA_draws = epRA_samples['ep_RA.{}.{}'.format(dna, constants[RBS])].sample(n_draws, replace=True).values
    ep_RA_mode = epRA_stats[epRA_stats['parameter']=='ep_RA.{}.{}'.format(dna, constants[RBS])]['mode'].values[0]
    for j, ind in enumerate(IND_muts):
        # Determine if ep_AI should be drawn. 
        if ind == 'Q294K':
            _kaki_samps = kaki_epAI_samples.sample(n_draws, replace=True)
            epAI_draws = _kaki_samps['ep_AI.{}'.format(ind)].values
            
            # Find the most-likely value for KaKi and ep_AI
            ka = kaki_epAI_stats[kaki_epAI_stats['parameter']=='Ka.{}'.format(ind)]['mode'].values[0]
            ki = kaki_epAI_stats[kaki_epAI_stats['parameter']=='Ki.{}'.format(ind)]['mode'].values[0]
            ep_ai = kaki_epAI_stats[kaki_epAI_stats['parameter']=='ep_AI.{}'.format(ind)]['mode'].values[0]
             
        else: 
            _kaki_samps = kaki_samples.sample(n_draws, replace=True)
            ka = kaki_stats[kaki_stats['parameter']=='Ka.{}'.format(ind)]['mode'].values[0]
            ki = kaki_stats[kaki_stats['parameter']=='Ki.{}'.format(ind)]['mode'].values[0] 
            epAI_draws = np.ones(n_draws) * constants['ep_AI']
            
            ep_ai = 4.5
            
        # Determine Ka/Ki
        ka_draws = _kaki_samps['Ka.{}'.format(ind)].values
        ki_draws = _kaki_samps['Ki.{}'.format(ind)].values
        
        # Compute the best fit
        fc = mut.thermo.SimpleRepression(R=constants[RBS], ep_r=ep_RA_mode,
                                        ka=ka, ki=ki, ep_ai=ep_ai, n_sites=constants['n_sites'],
                                        n_ns=constants['Nns'], effector_conc=c_range).fold_change()
        
        # Compute the credible region
        cred_region = np.zeros([2, len(c_range)])
        bohr_cred = np.zeros([2, len(c_range)])
        bohr_median = np.zeros(len(c_range))
        median = np.zeros(len(c_range))
        for k, c in enumerate(c_range):
            arch = mut.thermo.SimpleRepression(R=constants[RBS], ep_r=epRA_draws,       
                                             ka=ka_draws, ki=ki_draws,ep_ai=epAI_draws,
                                            n_sites=constants['n_sites'], n_ns=constants['Nns'],
                                            effector_conc=c)
            _fc = arch.fold_change()
            median[k] = np.median(_fc)
            cred_region[:, k] = mut.stats.compute_hpd(_fc, 0.95)
            bohr = arch.bohr_parameter()
            bohr_cred[:, k] = mut.stats.compute_hpd(bohr, 0.95)
            bohr_median[k] = np.median(bohr)
            
         # Assemble the dataframe. 
        _df = pd.DataFrame(np.array([c_range, fc, cred_region[0, :], cred_region[1, :],
                                        bohr_cred[0, :], bohr_cred[1, :], bohr_median]).T, 
                                       columns=['IPTGuM', 'mode', 'hpd_min', 'hpd_max', 
                                                'bohr_min', 'bohr_max', 'bohr_median'])
        _df['DNA_mutant'] = dna
        _df['IND_mutant'] = ind 
        dfs.append(_df)
            
df = pd.concat(dfs, sort=False)
df.to_csv('../../data/csv/Fig{}_{}_double_samples.csv'.format(FIG_NO, OPERATOR))
