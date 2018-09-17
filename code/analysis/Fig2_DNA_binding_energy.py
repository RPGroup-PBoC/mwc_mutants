# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import pystan
sys.path.insert(0, '../../')
import mut.bayes
import mut.stats
import imp
imp.reload(mut)
constants = mut.thermo.load_constants()
# Constants for file naming. 
FIG_NO = 2
OPERATOR = 'O2'

# %% Data loading & cleaning
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[(data['class'] == 'DNA') & (data['operator'] == 'O2') &
           (data['fold_change'] >=-0.2) & (data['fold_change'] <= 1.2)]
            
# Assign proper identifiers to each trace. 
DNA_idx = {m: i+1 for i, m in enumerate(data['mutant'].unique())}
REP_idx = {r:i+1 for i, r in enumerate(data['repressors'].unique())}
counter = 0
for m, i in DNA_idx.items(): 
    for r, idx in REP_idx.items():
        data.loc[(data['repressors']==r)& 
                     (data['mutant']==m), 'prof_idx'] = idx + counter
    counter += len(REP_idx)
data['prof_idx'] = data['prof_idx'].values.astype(int)

#%% Load and compile the inference model
model_code = mut.bayes.assemble_StanModelCode('../stan/DNA_binding_energy_induction.stan', '../stan/functions.stan')
model = pystan.StanModel(model_code=model_code)

# %% Assemble the data dictionary and perform inference.
data_dict = dict(J=len(data['prof_idx'].unique()), N=len(data), 
                idx=data['prof_idx'], R=data['repressors'],
                Nns=constants['Nns'], n_sites=constants['n_sites'],
                Ka=constants['Ka'], Ki=constants['Ki'], c=data['IPTGuM'],
                fc=data['fold_change'], ep_ai=constants['ep_AI'])
print('beginning sampling...')
samples = model.sampling(data_dict, iter=5000, chains=4, pars=['ep_RA', 'sigma'])
print('finished!')

# %% Clean sampling traces and format data frame
# Properly rename the parameters
new_names = {'ep_RA.{}'.format(i+1):'ep_RA.{}.{}'.format(DNA_data[DNA_data['prof_idx']==i+1]['mutant'].unique()[0],
                                                      int(DNA_data[DNA_data['prof_idx']==i+1]['repressors'].unique()[0])
                                                        ) for i in range(len(DNA_data['prof_idx'].unique()))}

# Format the dataframe. 
samples_df = mut.bayes.chains_to_dataframe(samples)
samples_df.rename(columns=new_names, inplace=True)

# Compute the statistics
samples_stats = mut.stats.compute_statistics(samples_df)
samples_df.to_csv('../../data/csv/Fig{}_{}_DNA_binding_energy_samples.csv'.format(FIG_NO, OPERATOR), index=False)
samples_stats.to_csv('../../data/csv/Fig{}_{}_DNA_binding_energy_stats.csv'.format(FIG_NO, OPERATOR), index=False)
print('Inference finished!')
