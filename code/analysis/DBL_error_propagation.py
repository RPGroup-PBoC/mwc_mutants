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

# Assign mutants identifiers. 
DNA_idx = {m: i+1 for i, m in enumerate(DNA_mutants)}
IND_idx = {m: i+1 for i, m in enumerate(IND_mutants)}
idx = {m: i+1 for i, m in enumerate(data['mutant'].unique())}
for i, dna in enumerate(DNA_mutants):
    for j, ind in enumerate(IND_mutants):
        data.loc[data['mutant'] == '{}-{}'.format(dna, ind), 'DNA_idx'] = DNA_idx[dna]
        data.loc[data['mutant'] == '{}-{}'.format(dna, ind), 'IND_idx'] = IND_idx[ind]
        data.loc[data['mutant'] == '{}-{}'.format(dna, ind), 'idx'] = idx['{}-{}'.format(dna, ind)]
data['DNA_idx'] = data['DNA_idx'].astype(int)
data['IND_idx'] = data['IND_idx'].astype(int)
data['idx'] = data['idx'].astype(int)

# Generate bounds.
ep_RA_upper = [stats_DNA[stats_DNA['parameter']=='ep_RA.{}'.format(d)]['hpd_max'].values[0] for d in DNA_idx.keys()]
ep_RA_lower = [stats_DNA[stats_DNA['parameter']=='ep_RA.{}'.format(d)]['hpd_min'].values[0] for d in DNA_idx.keys()]
ka_lower = [stats_IND[stats_IND['parameter']=='Ka.{}'.format(d)]['hpd_min'].values[0] for d in IND_idx.keys()]
ka_upper = [stats_IND[stats_IND['parameter']=='Ka.{}'.format(d)]['hpd_max'].values[0] for d in IND_idx.keys()]
ki_lower = [stats_IND[stats_IND['parameter']=='Ki.{}'.format(d)]['hpd_min'].values[0] for d in IND_idx.keys()]
ki_upper = [stats_IND[stats_IND['parameter']=='Ki.{}'.format(d)]['hpd_max'].values[0] for d in IND_idx.keys()]
print(ka_upper, ka_lower)

#%% Load and compile the stan model
model_code = mut.bayes.assemble_StanModelCode('../stan/DBL_error_propagation.stan', 
                                              '../stan/functions.stan')

print(model_code)
model = pystan.StanModel(model_code=model_code)
print('model compiled')

# %% Assemble the data dictionary and sample
data_dict = dict(J_DNA=len(DNA_mutants), J_IND=len(IND_mutants), J=len(data['mutant'].unique()),
                 N=len(data), DNA_idx=data['DNA_idx'].values, IND_idx=data['IND_idx'].values, 
                 idx=data['idx'].values, Nns=4.6E6, R=260, n_sites=2, ep_AI=4.5, c=data['IPTGuM'].values,
                 ep_RA_lower=ep_RA_lower, ep_RA_upper=ep_RA_upper, Ka_upper=ka_upper, Ka_lower=ka_lower,
                 Ki_upper=ki_upper, Ki_lower=ki_lower, fc=data['fold_change'].values)
print(stats_DNA, stats_IND)
print('sampling...')
samples = model.sampling(data_dict, iter=500, chains=4)
print('finished!')

# %%
samples_df = mut.bayes.chains_to_dataframe(samples)
new_names = {'ep_RA.{}'.format(i):'ep_RA.{}'.format(m) for m, i in DNA_idx.items()}
new_names.update({'Ka.{}'.format(i):'Ka.{}'.format(m) for m, i in DNA_idx.items()})
new_names.update({'Ki.{}'.format(i):'Ki.{}'.format(m) for m, i in DNA_idx.items()})
samples_df.rename(columns=new_names, inplace=True)
samples_df.to_csv('../../data/csv/DBL_binding_energy_samples.csv', index=False)
print('finished!')