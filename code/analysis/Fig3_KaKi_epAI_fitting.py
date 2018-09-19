# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import pystan
sys.path.insert(0, '../../')
import mut.bayes
import mut.stats
constants = mut.thermo.load_constants()

# Define experimental constants for inference. 
FIG_NO = 3
OPERATOR = 'O2'

# Data loading & cleaning
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[(data['class'] == 'IND') & (data['operator'] == 'O2') &
           (data['fold_change'] >=-0.2) & (data['fold_change'] <= 1.2)]


# Assign identifiers to data. 
IND_idx = {m:i+1 for i, m in enumerate(data['mutant'].unique())}
for m, idx in IND_idx.items():
    data.loc[data['mutant']==m, 'idx'] = idx
data['idx'] = data['idx'].astype(int)

# Assemble and load the inference model
model_code = mut.bayes.assemble_StanModelCode('../stan/KaKi_epAI_fitting.stan', '../stan/functions.stan')
model = pystan.StanModel(model_code=model_code)

# Assemble the data dictionary and sample the posterior.
data_dict = dict(J=len(IND_idx), N=len(data), idx=data['idx'], R=data['repressors'],
                Nns=constants['Nns'], c=data['IPTGuM'], ep_RA=constants[OPERATOR],
                n_sites=constants['n_sites'],fc=data['fold_change'])
samples = model.sampling(data_dict, iter=5000, chains=4)

# Clean the dataframe and rename variables. 
new_names = {'{}.{}'.format(n, i): '{}.{}'.format(n, m) for m, i in IND_idx.items() for n in ['Ka', 'Ki', 'ep_AI', 'sigma']}
samples_df = mut.bayes.chains_to_dataframe(samples)
samples_df.rename(columns=new_names, inplace=True)
samples_stats = mut.stats.compute_statistics(samples_df)
samples_df.to_csv('../../data/csv/Fig{}_{}_KaKi_epAI_samples.csv'.format(FIG_NO, OPERATOR), index=False)
samples_stats.to_csv('../../data/csv/Fig{}_{}_KaKi_epAI_stats.csv'.format(FIG_NO, OPERATOR), index=False)






