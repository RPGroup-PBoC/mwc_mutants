# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import mut.thermo
import mut.bayes
import mut.stats
import tqdm
constants = mut.thermo.load_constants()

# Load the prior predictive check data. 
prior_data = pd.read_csv('../../data/csv/epRA_prior_predictive_checks.csv')

# Load the stan model. 
model = mut.bayes.StanModel('../stan/DNA_binding_energy_induction.stan')

# Set up a dataframe to store the properties.
samples_dfs = []
sbc_dfs = []

# Definie the thinning constant for computing the rank statistic. 
thin = 5

# Iterate through each simulation
for g, d in tqdm.tqdm(prior_data.groupby('sim_idx')):
    
    # Determine the ground truth for each parameter.
    gt = {'ep_RA': d['ep_RA'].unique(),
         'sigma': d['sigma'].unique()}
    
    # Generate the data dictionary. 
    data_dict = {'J':1,
                'N': len(d),
                'idx': np.ones(len(d)).astype(int),
                'R': np.ones(len(d)) * constants['RBS1027'],
                'Nns': 4.6E6,
                'ep_ai': constants['ep_AI'],
                'n_sites': constants['n_sites'],
                'Ka': constants['Ka'],
                'Ki': constants['Ki'],
                'c': d['IPTGuM'],
                'fc': d['fc_draw']}
    
    # Sample the model
    _, samples = model.sample(data_dict=data_dict, iter=2000, control=dict(adapt_delta=0.99))
    samples.rename(columns={'ep_RA[1]': 'ep_RA', 'sigma[1]':'sigma'},
                  inplace=True)
    samples['sim_idx'] = g
    samples_dfs.append(_df)
    
    # Compute the properties for each parameter. 
    _sbc_dfs = []
    for p in ['ep_RA', 'sigma']:
        _df = pd.DataFrame([])
        z_score = (np.mean(samples[p]) - gt[p]) / np.std(samples[p])
        shrinkage = 1 - (np.var(samples[p]) / np.var(prior_data[p].unique()))
        _df['z_score'] = z_score
        _df['shrinkage'] = shrinkage
        _df['param'] = p 
        _df['rank'] = np.sum(samples[p].values[::thin] < gt[p])
        _df['ground_truth'] = gt[p]
        _sbc_dfs.append(_df)
        
    _sbc_dfs = pd.concat(_sbc_dfs)
    _sbc_dfs['sim_idx'] = g
    sbc_dfs.append(_sbc_dfs) 
sbc_df = pd.concat(sbc_dfs) 
    
sbc_df.to_csv('../../data/csv/epRA_sbc.csv', index=False)
