# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import pystan
sys.path.insert(0, '../../')
import mut.bayes
import mut.stats

# Set the constants.
N_ns = 4.6E6  # in base pairs
K_a = 139  # in µM
K_i = 0.53  # in µM
ep_ai = 4.5  # in k_BT.
n_sites = 2
# Load the data file.
data = pd.read_csv('../../data/csv/compiled_data.csv')
DNA = data[(data['class'] == 'DNA') | (data['class'] == 'WT')].copy()

# Include the identifier
idx_dict = {i: j for i, j in zip(
    DNA['mutant'].unique(), np.arange(1, 1 + len(DNA['mutant'].unique())))}
idx_key = {v: k for k, v in idx_dict.items()}

for m in DNA['mutant'].unique():
    DNA.loc[DNA['mutant'] == m, 'idx'] = idx_dict[m]

# Isolate the leakiness data.
leakiness = DNA[DNA['IPTGuM'] == 0].copy()

# Load the DNA binding energy stan model.
epR_model_code = mut.bayes.assemble_StanModelCode('../stan/hierarchical_epR_fit.stan',
                                                  '../stan/functions.stan')
epR_model = pystan.StanModel(model_code=epR_model_code)

# assemble the data dictionary.
data_dict = {'J': len(leakiness['mutant'].unique()), 'N': len(
    leakiness), 'trial': leakiness['idx'].values.astype(int),
    'R': leakiness['repressors'], 'c': leakiness['IPTGuM'], 'n_ns': N_ns, 'ka': K_a, 'ki': K_i,
    'ep_AI': ep_ai, 'n_sites': n_sites, 'fc': leakiness['fold_change']}

# Sample the posterior.
epR_chains = epR_model.sampling(data=data_dict, iter=10000, chains=4)
epR_df = mut.bayes.chains_to_dataframe(epR_chains)

# Rename the columns and save.
keys = epR_df.keys()
new_cols = {}
for i, k in enumerate(keys):
    if k != 'logp':
        param = k.split('.')
        new_cols[k] = '{}_{}'.format(param[0], idx_key[int(param[1])])

epR_df.rename(columns=new_cols, inplace=True)
epR_df.to_csv('../../data/mcmc/DNA_O2_epR_chains.csv', index=False)

# Compute and save the statistics.
epR_fit_statistics = mut.stats.compute_statistics(epR_df)
epR_fit_statistics.to_csv('../../data/mcmc/DNA_O2_epR_fit_statistics.csv')


# Load the stan model.
global_model_code = mut.bayes.assemble_StanModelCode(
    '../stan/hierarchical_global_fit.stan', '../stan/functions.stan')
global_model = pystan.StanModel(model_code=global_model_code)

# Define the data dictionary.
data_dict = {'J': len(DNA['mutant'].unique()),
             'N': len(DNA), 'trial': DNA['idx'].values.astype(int), 'R':DNA['repressors'],
             'n_ns':N_ns, 'c':DNA['IPTGuM'], 'ep_AI':ep_ai,
             'n_sites':n_sites, 'fc':DNA['fold_change'] }

# Sample the posterior.
global_fit_chains = global_model.sampling(data=data_dict, iter=10000, chains=4)
global_fit_df = mut.bayes.chains_to_dataframe(global_fit_chains)

# Rename the columns and save.
keys = global_fit_df.keys()
new_cols = {}
for i, k in enumerate(keys):
    if k != 'logp':
        param = k.split('.')
        new_cols[k] = '{}_{}'.format(param[0], idx_key[int(param[1])])

global_fit_df.rename(columns=new_cols, inplace=True)
global_fit_df.to_csv('../../data/mcmc/DNA_O2_global_fit_chains.csv', index=False)

# Compute and save the statistics.
global_fit_statistics = mut.stats.compute_statistics(global_fit_df)
global_fit_statistics.to_csv('../../data/mcmc/DNA_O2_global_fit_statistics.csv')
