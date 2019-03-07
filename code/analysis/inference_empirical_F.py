# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np 
import pandas as pd 
import mut.stats 
import mut.bayes
import mut.thermo
import tqdm
constants = mut.thermo.load_constants()

# Load the compiled data
data = pd.read_csv('../../data/csv/compiled_data.csv')
data.dropna(inplace=True)
data = data[(data['fold_change'] <= 1.2) & (data['fold_change'] >= -0.2)]

# Compute the reference bohr. 
ops = [constants[op] for op in data['operator']]
wt_bohr = mut.thermo.SimpleRepression(R=data['repressors'], ep_r=ops, 
                                       ka=constants['Ka'], ki=constants['Ki'],
                                       ep_ai=constants['ep_AI'], 
                                       effector_conc=data['IPTGuM']).bohr_parameter()
data['ref_bohr'] = wt_bohr

# Load the stan model. 
model = mut.bayes.StanModel('../stan/empirical_F.stan') #, force_compile=True)

# Assign unique identifiers. 
idx = data.groupby(['mutant', 'repressors', 'operator', 'IPTGuM']).ngroup() + 1 
data['idx'] = idx
data.sort_values('idx', inplace=True)

fc_stats = []
for g, d in tqdm.tqdm(data.groupby(['mutant', 'repressors', 'operator', 'IPTGuM'])):
    ref = d['ref_bohr'].unique()[0]
    # Assemble the data dictionary. 
    data_dict = {'N':len(d),
                  'ref_bohr':d['ref_bohr'].unique()[0],
                 'foldchange': d['fold_change']}
    fit, samples = model.sample(data_dict, iter=5000, control=dict(adapt_delta=0.99))
    
    dF_dfc = (samples['fc_mu'].median() - samples['fc_mu'].median()**2)**-1
    corr = dF_dfc * np.sign(ref) * ((samples['fc_sigma'].median() - samples['fc_sigma'].median()**2) - (np.exp(-ref)/(1 +    
                                                                                                    np.exp(-ref))**2)) 
    # Compute the stats. 
    _dbohr_stats = mut.stats.compute_statistics(samples, varnames=['delta_bohr', 'empirical_bohr', 'fc_mu', 'fc_sigma'], 
                                           logprob_name='lp__')    
    _dbohr_stats['correction'] = corr
    _dbohr_stats['mutant'] = g[0]
    _dbohr_stats['repressors'] = g[1]
    _dbohr_stats['operator'] = g[2]
    _dbohr_stats['IPTGuM'] = g[3]
    _dbohr_stats['class'] = d['class'].unique()[0]
    fc_stats.append(_dbohr_stats)
fc_stats = pd.concat(fc_stats)

fc_stats.to_csv('../../data/csv/empirical_F_statistics.csv', index=False)

