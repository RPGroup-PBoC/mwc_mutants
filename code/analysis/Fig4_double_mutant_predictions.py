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
RBS = 'RBS1027'

# %% Load the sampling traces. 
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[data['class']=='DBL']
epRA_samples = pd.read_csv('../../data/csv/DNA_binding_energy_samples.csv')
epRA_stats = pd.read_csv('../../data/csv/DNA_binding_energy_summary.csv')
kaki_samples = pd.read_csv('../../data/csv/KaKi_only_samples.csv')
kaki_stats = pd.read_csv('../../data/csv/KaKi_only_summary.csv')
kaki_epAI_samples = pd.read_csv('../../data/csv/KaKi_epAI_samples.csv')
kaki_epAI_stats = pd.read_csv('../../data/csv/KaKi_epAI_summary.csv')

# Define the sampling constants. 
n_draws = int(1E5)

# Get the mutants
muts = data['mutant'].unique()
IND_muts = ['Q294K', 'F164T', 'Q294V']
DNA_muts = ['Y20I', 'Q21M', 'Q21A']

# Loop through mutants and isolate parameters.
dfs = []
c_range = np.logspace(-2, 4, 200)
for i, dna in enumerate(tqdm.tqdm(DNA_muts)):
    epRA_draws = epRA_samples[(epRA_samples['mutant']==dna) & (epRA_samples['repressors']==260) & 
                              (epRA_samples['operator']=='O2')]['ep_RA'].sample(n_draws, replace=True).values
    ep_RA_mode = epRA_stats[(epRA_stats['parameter']=='ep_RA') & (epRA_stats['mutant']==dna) &
                            (epRA_stats['operator']=='O2') & (epRA_stats['repressors']==260)]['mode'].values[0]
    for j, ind in enumerate(IND_muts):
        # Determine if ep_AI should be drawn. 
        if ind == 'Q294K':
            _kaki_samps = kaki_epAI_samples[(kaki_epAI_samples['mutant']==ind) & (kaki_epAI_samples['operator']=='O2')].sample(n_draws, replace=True)
            _kaki_stats = kaki_epAI_stats[(kaki_epAI_stats['mutant']==ind) & (kaki_epAI_stats['operator']=='O2')] 
            ep_ai = _kaki_stats[_kaki_stats['parameter']=='ep_AI']['mode'].values[0]
        else: 
            _kaki_samps = kaki_samples[(kaki_samples['mutant']==ind) & (kaki_samples['operator']=='O2')].sample(n_draws, replace=True)
            _kaki_samps['ep_AI'] = constants['ep_AI']
            _kaki_stats = kaki_stats[(kaki_stats['mutant']==ind) & (kaki_stats['operator']=='O2')] 
            ep_ai = constants['ep_AI']
        
        # Find the most-likely value for KaKi and ep_AI
        ka = _kaki_stats[_kaki_stats['parameter']=='Ka']['mode'].values[0]
        ki = _kaki_stats[_kaki_stats['parameter']=='Ki']['mode'].values[0]  
       
        # Determine Ka/Ki
        ka_draws = _kaki_samps['Ka'].values
        ki_draws = _kaki_samps['Ki'].values
        epAI_draws = _kaki_samps['ep_AI'].values
        
        # Compute the best fit
        fc = mut.thermo.SimpleRepression(R=260, ep_r=ep_RA_mode,
                                        ka=ka, ki=ki, ep_ai=ep_ai, n_sites=constants['n_sites'],
                                        n_ns=constants['Nns'], effector_conc=c_range).fold_change()
        
        # Compute the credible region
        cred_region = np.zeros([2, len(c_range)])
        bohr_cred = np.zeros([2, len(c_range)])
        bohr_median = np.zeros(len(c_range))
        median = np.zeros(len(c_range))
        for k, c in enumerate(c_range):
            arch = mut.thermo.SimpleRepression(R=260, ep_r=epRA_draws,       
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
df.to_csv('../../data/csv/DBL_mutant_predictions.csv')
