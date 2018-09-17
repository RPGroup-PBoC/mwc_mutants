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
FIG_NO = 1
OPERATOR = 'O2'

# Load and trim the data
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[(data['class'] == 'DNA') & (data['operator'] == 'O2') &
           (data['fold_change'] >=-0.2) & (data['fold_change'] <= 1.2)]

# Assign identifiers to each mutant. 
DNA_idx = {m:i+1 for i, m in enumerate(data['mutant'].unique())}
for m, i in DNA_idx.items():
    data.loc[data['mutant']==m, 'idx'] = DNA_idx[m]
data['idx'] = data['idx'].astype(int)

# Assemble and load the model code
model_code = mut.bayes.assemble_StanModelCode('../stan/DNA_binding_energy_induction.stan', '../stan/functions.stan')
model = pystan.StanModel(model_code=model_code)

# Perform the inference
data_dict = dict(J=len(data['idx'].unique()), N=len(data), 
                idx=data['idx'], R=data['repressors'],
                Nns=constants['Nns'], n_sites=constants['n_sites'],
                Ka=constants['Ka'], Ki=constants['Ki'], c=data['IPTGuM'],
                fc=data['fold_change'], ep_ai=constants['ep_AI'])
samples = model.sampling(data_dict, iter=5000, chains=4)
print(samples)
# Properly rename the parameters
new_names = {'{}.{}'.format(n, i):'{}.{}'.format(n, m) for m, i in DNA_idx.items() for n in ['ep_RA', 'sigma']}

# Format the dataframe. 
samples_df = mut.bayes.chains_to_dataframe(samples)
samples_df.rename(columns=new_names, inplace=True)

# Compute the statistics
samples_stats = mut.stats.compute_statistics(samples_df)
samples_df.to_csv('../../data/csv/FigS{}_{}_DNA_binding_energy_global_samples.csv'.format(FIG_NO, OPERATOR), index=False)
samples_stats.to_csv('../../data/csv/FigS2{}_{}_DNA_binding_energy_global_stats.csv'.format(FIG_NO, OPERATOR), index=False)
print('Inference finished!')
