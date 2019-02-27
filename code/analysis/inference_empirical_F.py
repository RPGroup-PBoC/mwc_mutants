# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np 
import pandas as pd 
import mut.stats 
import mut.bayes
import mut.thermo
constants = mut.thermo.load_constants()

# Load the compiled data
data = pd.read_csv('../../data/csv/compiled_data.csv')
data.dropna(inplace=True)

# Compute the reference bohr. 
ops = [constants[op] for op in data['operator']]
wt_bohr = -mut.thermo.SimpleRepression(R=data['repressors'], ep_r=ops, 
                                       ka=constants['Ka'], ki=constants['Ki'],
                                       ep_ai=constants['ep_AI'], 
                                       effector_conc=data['IPTGuM']).bohr_parameter()
data['ref_bohr'] = wt_bohr

# Load the stan model. 
model = mut.bayes.StanModel('../stan/empirical_F.stan', force_compile=True)

# Assign unique identifiers. 
idx = data.groupby(['mutant', 'IPTGuM', 'repressors', 'operator']).ngroup() + 1 
data['idx'] = idx

# Assemble the data dictionary. 
data_dict = {'N':len(data),
             'J':data['idx'].max(),
             'idx':data['idx'], 
             'foldchange': data['fold_change'],
             'bohr_ref':data.groupby(['idx'])['ref_bohr'].mean().values}

# Sample
fit, samples = model.sample(data_dict, iter=1000) 


# Compute the statistics and add identifiers. 
fc_vars = [f'fc_mu[{i}]' for i in idx.values]
bohr_vars = [f'empirical_bohr[{i}]' for i in idx.values]
delta_bohr_vars = [f'delta_bohr[{i}]' for i in idx.values]
fc_stats = mut.stats.compute_statistics(samples, varnames=fc_vars, 
                                        logprob_name='lp__')
bohr_stats = mut.stats.compute_statistics(samples, varnames=bohr_vars,
                                          logprob_name='lp__')
dbohr_stats = mut.stats.compute_statistics(samples, varnames=bohr_vars,
                                          logprob_name='lp__')

for d in [fc_stats, bohr_stats, dbohr_stats]:
    d['mutant'] = data['mutant'].values
    d['IPTGuM'] = data['IPTGuM'].values
    d['repressors'] = data['repressors'].values
    d['operator'] = data['operator'].values
    d['class'] = data['class'].values

fc_stats.drop(columns=['parameter'], inplace=True)
fc_stats.to_csv('../../data/csv/empirical_F_statistics.csv', index=False)

# Compute the posterior predictive checks
rep_vars = [f'y_rep[{i+1}]' for i in range(len(data))]
rep_stats = mut.stats.compute_statistics(samples, varnames=rep_vars, 
                                         logprob_name='lp__')
rep_stats.to_csv('../../data/csv/inferred_fc_ppc.csv', index=False)

