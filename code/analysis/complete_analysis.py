# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import pandas as pd
import pystan
sys.path.insert(0, '../../')
import mut.bayes
import mut.stats

# %%
# Filter the data and exclude nonsense.
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[(data['fold_change'] >= -0.1) & (data['fold_change'] <= 1.2)]
data.dropna(inplace=True, axis=0)

# Separate into the relevant classes.
DNA = data[data['class']=='DNA']
IND = data[data['class']=='IND']
DBL = data[data['class']=='DBL']

# Determine the number of unique mutants in each class.
J_DNA = len(DNA['mutant'].unique())
J_IND = len(IND['mutant'].unique())
J_DBL = len(DBL['mutant'].unique())

# Get the number of measurements for each class
N_DNA = len(DNA)
N_IND = len(IND)
N_DBL = len(DBL)

# Define the identifiers and save as a key. 
DNA_key = {m:int(i+1) for i, m in enumerate(DNA['mutant'].unique())}
IND_key = {m:int(i+1) for i, m in enumerate(IND['mutant'].unique())}
DBL_key = {m:int(i+1) for i, m in enumerate(DBL['mutant'].unique())}
for m in DNA['mutant'].unique():
    DNA.loc[DNA['mutant'] == m, 'idx'] = DNA_key[m]
for m in IND['mutant'].unique():
    IND.loc[IND['mutant'] == m, 'idx'] = IND_key[m]
for m in DBL['mutant'].unique():
    DBL.loc[DBL['mutant'] == m, 'idx'] = DBL_key[m]

# Isolate the identifier vectors.
idx_DNA = DNA['idx'].values.astype(int)
idx_IND = IND['idx'].values.astype(int)
idx_DBL = DBL['idx'].values.astype(int)

# Define repressor copy number vectors. 
R_DNA = DNA['repressors']
R_IND = IND['repressors']
R_DBL = DBL['repressors']

# Define the fold-change measurements.
fc_DNA = DNA['fold_change']
fc_IND = IND['fold_change']
fc_DBL = DBL['fold_change']

# Define the concentration vectors.
c_DNA = DNA['IPTGuM']
c_IND = IND['IPTGuM']
c_DBL = DBL['IPTGuM']

# Define the thermodynamic model constants.
ep_ai = 4.5 # in units of kBT
Nns = 4.6E6 # Length of E. coli genome
n_sites = 2
wt_epR = -13.9 # in units of kBT
wt_ka = 139 # in units of µM
wt_ki = 0.53 # in units of µM


# Assemble the complete data dictionary.
data_dict = {'J_DNA':J_DNA, 'J_IND':J_IND, 'J_DBL':J_DBL, 
             'N_DNA':N_DNA, 'N_IND':N_IND, 'N_DBL':N_DBL,
             'idx_DNA':idx_DNA, 'idx_IND':idx_IND, 'idx_DBL':idx_DBL, 
             'R_DNA':R_DNA, 'R_IND':R_IND, 'R_DBL':R_DBL,
             'wt_ka':wt_ka, 'wt_ki':wt_ki, 'wt_epR':wt_epR,
             'ep_ai':ep_ai, 'Nns':Nns, 'n_sites':int(n_sites),
             'c_DNA':c_DNA, 'c_IND':c_IND, 'c_DBL':c_DBL,
             'fc_DNA':fc_DNA, 'fc_IND':fc_IND, 'fc_DBL':fc_DBL}

# %% Compile and execute the stan model.
model = pystan.StanModel('../stan/complete_parameter_estimation.stan')

# %%
samples = model.sampling(data=data_dict, iter=5000, chains=2)