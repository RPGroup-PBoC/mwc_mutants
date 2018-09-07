# -*- coding: utf-8 -*-
#%%
import sys
import numpy as np
import pandas as pd
import pystan
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.bayes
import mut.stats

# %% Data loading & cleaning
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[(data['class'] == 'DNA') & (data['fold_change'] >= 0) &\
            (data['operator'] == 'O2') & (data['IPTGuM'] == 0)]

# Assign identifiers to unique mutants.
DNA_idx = {m: i+1 for i, m in enumerate(data['mutant'].unique())}
for m, idx in DNA_idx.items():
    data.loc[data['mutant'] == m, 'idx'] = DNA_idx[m]
data['idx'] = data['idx'].astype('int')

#%% Load and compile the inference model
model_code = mut.bayes.assemble_StanModelCode('../stan/DNA_binding_energy.stan', '../stan/functions.stan')
print(model_code)
model = pystan.StanModel(model_code=model_code)
print('Compilation finished!')

# %% Assemble the data dictionary and perform inference.
data_dict = dict(J=len(DNA_idx), N=len(data), idx=data['idx'], R=data['repressors'], Nns=4.6E6,
                 ep_ai=4.5, n_sites=2, fc=data['fold_change'])
print('beginning sampling...')
samples = model.sampling(data_dict, iter=10000, chains=4, pars=['ep_RA', 'sigma'])
print('finished!')

# %% Clean sampling traces and format data frame
samples_df = mut.bayes.chains_to_dataframe(samples)
new_names = {'ep_RA.{}'.format(i):'ep_RA.{}'.format(m) for m, i in DNA_idx.items()}
samples_df.rename(columns=new_names, inplace=True)
samples_df.to_csv('../../data/csv/DNA_binding_energy_samples.csv', index=False)

