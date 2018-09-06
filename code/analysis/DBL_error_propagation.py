#%%
import sys
import numpy as np 
import pandas as pd
import pystan
sys.path.insert(0, '../../')
import mut.bayes
import mut.stats

# Load the experimental data.
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[(data['class']=='DBL') & (data['fold_change'] >= 0) &\
            (data['fold_change'] <= 1.2) & (data['operator'] == 'O2')].copy()

# Define a dictionary of the measured doulbe mutants separated by domain.
DNA_mutants = list(set([m.split('-')[0] for m in data['mutant'].unique()]))
IND_mutants = list(set([m.split('-')[1] for m in data['mutant'].unique()]))

# Load the parameter inference data.
samples_DNA = pd.read_csv('../../data/csv/20180902_DNA_binding_energy_samples.csv')
stats_DNA = mut.stats.compute_statistics(samples_DNA)
samples_IND = pd.read_csv('../../data/csv/20180902_KaKi_samples.csv')
stats_IND = mut.stats.compute_statistics(samples_IND)

#%% Load and compile the stan model
model_code = mut.bayes.assemble_StanModelCode('../stan/DBL_error_propagation.stan', 
                                              '../stan/functions.stan')
