# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import pystan
import glob
sys.path.insert(0, '../../')
import mut.bayes
import mut.stats

# Set the constants.
N_ns = 4.6E6  # in base pairs
ep_ai = 4.5  # in k_BT.
n_sites = 2

# Load the data
data = pd.read_csv('../../data/csv/compiled_data.csv')
DBL = data[(data['class'] == 'DBL') | (data['class'] == 'WT')].copy()

# Add the mutant identifiers.
idx_dict = {i: j for i, j in zip(
    DBL['mutant'].unique(), np.arange(1, 1 + len(DBL['mutant'].unique())))}
idx_key = {v: k for k, v in idx_dict.items()}

for m in DBL['mutant'].unique():
    DBL.loc[DBL['mutant'] == m, 'idx'] = idx_dict[m]

# Load the stan modeled for the global fit.
global_fit_model_code = mut.bayes.assemble_StanModelCode(
    '../stan/hierarchical_global_fit.stan', '../stan/functions.stan')
global_fit_model = pystan.StanModel(model_code=global_fit_model_code)

# Assemble the dataframe and sample the distribution.
data_dict = {'J': len(DBL['mutant'].unique()), 'N': len(
    DBL), 'trial': DBL['idx'].values.astype(int), 'R':DBL['repressors'],
    'n_ns':N_ns, 'ep_AI':ep_ai, 'n_sites':n_sites, 'c':DBL['IPTGuM'],
    'fc':DBL['fold_change']}
global_fit_chains = global_fit_model.sampling(data=data_dict, iter=10000)
global_fit_df = mut.bayes.chains_to_dataframe(global_fit_chains)

# Rename the columns and save.
keys = global_fit_df.keys()
new_cols = {}
for i, k in enumerate(keys):
    if k != 'logp':
        param = k.split('.')
        new_cols[k] = '{}_{}'.format(param[0], idx_key[int(param[1])])

global_fit_df.rename(columns=new_cols, inplace=True)
global_fit_df.to_csv('../../data/mcmc/DBL_O2_global_fit_chains.csv', index=False)

global_fit_stats = mut.stats.compute_statistics(global_fit_df)
global_fit_stats.to_csv('../../data/mcmc/DBL_O2_global_fit_statistics.csv', index=False)