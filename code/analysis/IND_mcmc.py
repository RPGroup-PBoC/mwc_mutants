# %%
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
ep_R = -13.9  # in k_BT
ep_ai = 4.5  # in k_BT.
n_sites = 2

# Load the data and format as desired.
data = pd.read_csv('../../data/csv/compiled_data.csv')
IND = data[(data['class'] == 'IND') | (data['class'] == 'WT')].copy()
IND = IND[IND['fold_change'] >= 0]

# Include the identifier
idx_dict = {i: j for i, j in zip(
    IND['mutant'].unique(), np.arange(1, 1 + len(IND['mutant'].unique())))}
idx_key = {v: k for k, v in idx_dict.items()}

for m in IND['mutant'].unique():
    IND.loc[IND['mutant'] == m, 'idx'] = idx_dict[m]
# Load the KaKi analysis stan model
KaKi_model_code = mut.bayes.assemble_StanModelCode('../stan/hierarchical_kaki_fit.stan',
                                                   '../stan/functions.stan')
KaKi_model = pystan.StanModel(model_code=KaKi_model_code)

#%% assemble the data dictionary.
data_dict = {'J': len(IND['mutant'].unique()), 'N': len(
    IND), 'trial': IND['idx'].values.astype(int),
    'R': IND['repressors'], 'c': IND['IPTGuM'], 'n_ns': N_ns, 'ep_R': ep_R,
    'ep_AI': ep_ai, 'n_sites': n_sites, 'fc': IND['fold_change']}

KaKi_chains = KaKi_model.sampling(
    data=data_dict, iter=10000, chains=4)
KaKi_df = mut.bayes.chains_to_dataframe(KaKi_chains)

# Rename the columns and save.
keys = KaKi_df.keys()
new_cols = {}
for i, k in enumerate(keys):
    if k != 'logp':
        param = k.split('.')
        new_cols[k] = '{}_{}'.format(param[0], idx_key[int(param[1])])

KaKi_df.rename(columns=new_cols, inplace=True)
KaKi_df.to_csv('../../data/mcmc/IND_O2_KaKi_fit_chains.csv', index=False)

stats = mut.stats.compute_statistics(KaKi_df)
stats.to_csv('../../data/mcmc/IND_O2_KaKi_fit_statistics.csv', index=False)


#%% Perform a global fit including the binding energy.
global_fit_model_code = mut.bayes.assemble_StanModelCode(
    '../stan/hierarchical_global_fit.stan', '../stan/functions.stan')
global_fit_model = pystan.StanModel(model_code=global_fit_model_code)

# Assemble the data set for the global fit model.
data_dict = {'J': len(IND['mutant'].unique()), 'N': len(
    IND), 'trial': IND['idx'].values.astype(int), 'R': IND['repressors'],
    'n_ns': N_ns, 'ep_AI': ep_ai, 'n_sites': n_sites, 'c': IND['IPTGuM'],
    'fc': IND['fold_change']}

# Sample and convert to a dtaframe.
global_fit_chains = global_fit_model.sampling(
    data=data_dict, chains=48, iter=50000, thin=50)
global_fit_df = mut.bayes.chains_to_dataframe(global_fit_chains)

# Rename the columns and save.
keys = KaKi_df.keys()
new_cols = {}
for i, k in enumerate(keys):
    if k != 'logp':
        param = k.split('.')
        new_cols[k] = '{}_{}'.format(param[0], idx_key[int(param[1])])
global_df.rename(coumns=new_cols, inplace=True)
global_fit_df.to_csv('../../data/mcmc/IND_O2_global_fit_chains.csv', index=False)

global_fit_stats = mut.stats.compute_statistics(global_fit_df)
global_fit_stats.to_csv('../../data/mcmc/IND_O2_global_fit_statistics.csv', index=False)

