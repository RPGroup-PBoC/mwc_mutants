# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np 
import pandas as pd 
import mut.stats 
import mut.bayes

# Load the compiled data
data = pd.read_csv('../../data/csv/compiled_data.csv')
data.dropna(inplace=True)

# Load the stan model. 
model = mut.bayes.StanModel('../stan/empirical_F.stan')

# Assign unique identifiers. 
idx = data.groupby(['mutant', 'IPTGuM', 'repressors', 'operator']).ngroup() + 1 
data['idx'] = idx

# Assemble the data dictionary. 
data_dict = {'N':len(data),
             'J':data['idx'].max(),
             'idx':data['idx'], 
             'foldchange': data['fold_change']}

# Sample
fit, samples = model.sample(data_dict, iter=10000,
                            control=dict(adapt_delta=0.999))


# Compute the statistics and add identifiers. 
fc_vars = [f'fc_mu[{i}]' for i in idx.values]
bohr_vars = [f'empirical_bohr[{i}]' for i in idx.values]
fc_stats = mut.stats.compute_statistics(samples, varnames=fc_vars, logprob_name='lp__')
fc_stats.rename(columns={'mean': 'fc_mean', 'median':'fc_median', 'mode':'fc_mode', 'hpd_min':'fc_min', 
                         'hpd_max':'fc_max'}, inplace=True)
bohr_stats = mut.stats.compute_statistics(samples, varnames=bohr_vars, logprob_name='lp__')
bohr_stats.rename(columns={'mean': 'bohr_mean', 'median':'bohr_median', 'mode':'bohr_mode', 'hpd_min':'bohr_min', 
                           'hpd_max':'bohr_max'}, inplace=True)

# Add Identifiying information. 
fc_stats['bohr_median'] = bohr_stats['bohr_median']
fc_stats['bohr_mean'] = bohr_stats['bohr_mean']
fc_stats['bohr_mode']=bohr_stats['bohr_mode']
fc_stats['bohr_min']=bohr_stats['bohr_min']
fc_stats['bohr_max']=bohr_stats['bohr_max']
fc_stats['mutant'] = data['mutant'].values
fc_stats['IPTGuM'] = data['IPTGuM'].values
fc_stats['repressors'] = data['repressors'].values
fc_stats['operator'] = data['operator'].values
fc_stats['class'] = data['class'].values
fc_stats.drop(columns=['parameter'], inplace=True)
fc_stats.to_csv('../../data/csv/empirical_F_statistics.csv', index=False)
