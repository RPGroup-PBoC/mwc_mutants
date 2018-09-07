# -*- coding:utf-8 -*-
#%%
import sys
import numpy as np 
import pandas as pd 
import pystan
sys.path.insert(0, '../../')
import mut.stats
import mut.bayes

# Load and prepare the data
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[(data['class']=='IND') & (data['operator']=='O2') &
            (data['fold_change'] >= 0)]
IND_idx = {m:i+1 for i, m in enumerate(data['mutant'].unique())}
for m, idx in IND_idx.items():
    data.loc[data['mutant']==m, 'idx'] = idx
data['idx'] = data['idx'].astype(int)

# %%  Load and compile the stan model. 
model_code = mut.bayes.assemble_StanModelCode('../stan/KaKi_fitting.stan',
                                               '../stan/functions.stan')
model = pystan.StanModel(model_code=model_code)
print(model_code)


# %% Assemble the data dictionary and sample.
data_dict = dict(J=len(IND_idx), N=len(data), idx=data['idx'], R=data['repressors'],
                 Nns=4.6E6, c=data['IPTGuM'], ep_RA=-13.9, ep_AI=4.5, n_sites=2,
                 fc=data['fold_change'])
samples = model.sampling(data_dict, iter=5000, chains=4, pars=['Ka', 'Ki', 'sigma'])

# %%
samples_df = mut.bayes.chains_to_dataframe(samples)
new_names = {'Ka.{}'.format(i):'Ka.{}'.format(m) for m, i in IND_idx.items()}
new_names.update({'Ki.{}'.format(i):'Ki.{}'.format(m) for m, i in IND_idx.items()})
samples_df.rename(columns=new_names, inplace=True)
samples_df.to_csv('../../data/csv/20180902_KaKi_samples.csv', index=False)

# %%