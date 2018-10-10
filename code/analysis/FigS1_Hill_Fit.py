# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import pystan
sys.path.insert(0, '../../')
import mut.stats
import mut.bayes

# Load the data 
data = pd.read_csv('../../data/csv/compiled_data.csv')

# Isolate the wt and DNA binding mutants. 
data = data[((data['class']=='WT') | (data['class']=='DNA')) & (data['operator']=='O2')].copy()

# Add run identifiers
idx = {}
counter = 1
for m, i in mut_idx.items():
    for r, z in rep_idx.items():
        data.loc[(data['mutant']==m) & (data['repressors']==r), 'idx'] = counter
        idx[f'{m.upper()}.{int(r)}'] = counter
        counter += 1    
data['idx'] = data['idx'].astype(int)

# Compile the inference model. 
model = pystan.StanModel('../stan/Hill_fitting.stan')

# Assemble the data dictionary and sample. 
data_dict = dict(J=len(idx), N=len(data), idx=data['idx'],
                c=data['IPTGuM'], foldchange=data['fold_change'])
samples = model.sampling(data_dict, iter=5000, chains=4)

# Generate dataframe, rename parameters, and save. 
samples_df = mut.bayes.chains_to_dataframe(samples)
new_names = {f'{p}.{i}':f'{p}.{m.upper()}' for m, i in idx.items() for p in ['a', 'b', 'K', 'n', 'sigma']}
samples_df.rename(columns=new_names, inplace=True)

# Compute the statistics. 
stats = mut.stats.compute_statistics(samples_df)

# Save to disk. 
samples_df.to_csv('../../data/csv/FigS1_Hill_samples.csv', index=False)
stats.to_csv('../../data/csv/FigS1_Hill_stats.csv', index=False)

idx
