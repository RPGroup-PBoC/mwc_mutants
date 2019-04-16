# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import mut.thermo
import mut.bayes
import mut.stats
import joblib
import tqdm
constants = mut.thermo.load_constants()

# Load the prior predictive check data. 
prior_data = pd.read_csv('../../data/csv/IND_prior_predictive_checks.csv')


# Load the stan model. 
KaKi_model = mut.bayes.StanModel('../stan/KaKi_fitting.stan')
KaKi_epAI_model = mut.bayes.StanModel('../stan/KaKi_epAI_fitting.stan') 
_model = {'KaKi_only': KaKi_model, 'KaKi_epAI':KaKi_epAI_model}

# Set up a dataframe to store the properties.
samples_dfs = []
sbc_dfs = []

# Definie the thinning constant for computing the rank statistic. 
thin = 5

# Iterate through each simulation
for g, d in tqdm.tqdm(prior_data.groupby(['model', 'draw'])):
    # Generate the data dictionary. 
    data_dict = {'J':1,
                'N': len(d),
                'idx': np.ones(len(d)).astype(int),
                'ep_RA': constants['O2'],
                'R': np.ones(len(d)) * constants['RBS1027'],
                'Nns': 4.6E6,
                'n_sites': constants['n_sites'],
                'c': d['IPTGuM'],
                'fc': d['fc_draw']}
    # Define the columns for renaming
    columns={'Ka[1]': 'Ka', 'sigma[1]':'sigma', 'Ki[1]':'Ki'}

    # Determine the ground truth for each parameter.
    if g[0] == 'KaKi_only':
        gt = {'Ka': d['ka'].unique(),
              'Ki': d['ki'].unique(),
              'sigma': d['sigma'].unique()}
        data_dict['ep_AI'] = constants['ep_AI']
        pars = ['Ka', 'Ki', 'sigma']
    else:
        gt = {'Ka': d['ka'].unique(),
              'Ki': d['ki'].unique(),
              'ep_AI': d['ep_ai'].unique(),
              'sigma': d['sigma'].unique()}

        pars = ['Ka', 'Ki', 'ep_AI', 'sigma']
        columns['ep_AI[1]'] = 'ep_AI'

        
    # Sample the model
    model = _model[g[0]]
    _, samples = model.sample(data_dict=data_dict, iter=2000, 
                control=dict(adapt_delta=0.99), n_jobs=1)

    samples.rename(columns,  inplace=True)
    samples['sim_idx'] = g[0]
    samples['model'] = g[1]
    samples_dfs.append(samples)
    
    # Compute the properties for each parameter. 
    _sbc_dfs = []
    for p in pars:
        _df = pd.DataFrame([])
        z_score = (np.mean(samples[p]) - gt[p]) / np.std(samples[p])
        shrinkage = 1 - (np.var(samples[p]) / np.var(prior_data[p].unique()))
        _df['z_score'] = z_score
        _df['shrinkage'] = shrinkage
        _df['param'] = p 
        _df['rank'] = np.sum(samples[p].values[::thin] < gt[p])
        _df['rank_ndraws'] = len(samples[p].values[::thin])
        _df['post_median'] = np.mean(samples[p])
        _df['post_mean'] = np.median(samples[p])
        _df['post_mode'] = samples.iloc[np.argmax(samples['lp__'].values)][p]
        _df['model'] = g[0]
        _df['ground_truth'] = gt[p]
        _sbc_dfs.append(_df)
        
    _sbc_dfs = pd.concat(_sbc_dfs)
    _sbc_dfs['sim_idx'] = g[1]
    sbc_dfs.append(_sbc_dfs) 
sbc_df = pd.concat(sbc_dfs) 
    
sbc_df.to_csv('../../data/csv/IND_sbc.csv', index=False)
samples_df = pd.concat(samples_dfs)
samples_df.to_csv('../../data/csv/IND_sbc_samples.csv', index=False)
