# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import pystan
sys.path.insert(0, '../../')
import mut.bayes
import mut.stats
constants = mut.thermo.load_constants()

# Constants for file naming. 
FIG_NO = 2
OPERATOR = 'O2'
RBS = 'RBS1027'

# Load the data
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[(data['class']=='DBL') & (data['fold_change'] >=-0.2) & 
           (data['fold_change'] <= 1.2)]

# Assign identifiers. 
idx = {m:i+1 for i, m in enumerate(data['mutant'].unique())}
for m, i in idx.items():
    data.loc[data['mutant']==m, 'idx'] = i
data['idx'] = data['idx'].astype(int)

# Assemble and compile the stan model. 
model_code = mut.bayes.assemble_StanModelCode('../stan/DBL_parameter_estimation.stan', 
                                              '../stan/functions.stan')
model =  pystan.StanModel(model_code = model_code)

# Assemble the data dictionary and sample. 
data_dict = dict(J=len(idx), N=len(data), idx=data['idx'], R=data['repressors'],
                Nns=constants['Nns'], n_sites=constants['n_sites'], 
                ep_AI=constants['ep_AI'], c=data['IPTGuM'],
                fc=data['fold_change'])
samples = model.sampling(data_dict, iter=7000, chains=4, pars=['Ka', 'Ki', 'ep_RA', 'sigma'])

# Convert to dataframe and clean up
samples_df = mut.bayes.chains_to_dataframe(samples)
new_names = {'{}.{}'.format(n, i):'{}.{}'.format(n, m) for m, i in idx.items() for n in ['ep_RA', 'Ka', 'Ki', 'sigma']}
samples_df.rename(columns=new_names, inplace=True)
samples_stats = mut.stats.compute_statistics(samples_df)
samples_df.to_csv('../../data/csv/Fig{}_{}_DBL_samples.csv'.format(FIG_NO, OPERATOR))
samples_stats.to_csv('../../data/csv/Fig{}_{}_DBL_stats.csv'.format(FIG_NO, OPERATOR))
