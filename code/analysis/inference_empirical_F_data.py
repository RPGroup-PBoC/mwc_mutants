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
data = data[data['mutant']!='wt']
ind_data = pd.read_csv('../../data/csv/RazoMejia2018_data.csv')
ind_data['repressors'] *= 2
ind_data.rename(columns={'fold_change_A':'fold_change',
                        'IPTG_uM':'IPTGuM'}, inplace=True)


# Load the stan model. 
model = mut.bayes.StanModel('../stan/empirical_F_data.stan') #, force_compile=True)


# Compute the reference bohr. 
ops = [constants[op] for op in data['operator']]
wt_bohr = -mut.thermo.SimpleRepression(R=data['repressors'], ep_r=ops, 
                                       ka=constants['Ka'], ki=constants['Ki'],
                                       ep_ai=constants['ep_AI'], 
                                       effector_conc=data['IPTGuM']).bohr_parameter()
data['ref_bohr'] = wt_bohr
# Assign unique identifiers. 
idx = data.groupby(['mutant', 'repressors', 'operator', 'IPTGuM']).ngroup() + 1 
data['idx'] = idx
data.sort_values('idx', inplace=True)

fc_stats = []
for g, d in tqdm.tqdm(data.groupby(['mutant', 'repressors', 'operator', 'IPTGuM'])):
    _wt = ind_data[(ind_data['repressors']==g[1]) &
                 (ind_data['operator']==g[2]) & 
                  (ind_data['IPTGuM']==g[3])]
    if len(_wt) > 0:
        # Assemble the data dictionary. 
        data_dict = {'N_mut':len(d),
                     'N_wt':len(_wt),
                     'fc_wt': _wt['fold_change'],
                     'ref_bohr': d['ref_bohr'].unique()[0],
                     'fc_mut': d['fold_change']}
        fit, sample = model.sample(data_dict) #, iter=2000, control=dict(adapt_delta=0.99))
        sample['bohr_wt'] = -np.log(sample['fc_mu_wt']**-1 -  1)
        sample['bohr_mut'] = -np.log(sample['fc_mu_mut']**-1 - 1)
        sample['delta_bohr'] = sample['bohr_wt'] - sample['bohr_mut']
        _dbohr_stats = mut.stats.compute_statistics(sample, logprob_name='lp__')
        _dbohr_stats['mutant'] = g[0]
        _dbohr_stats['repressors'] = g[1]
        _dbohr_stats['operator'] = g[2]
        _dbohr_stats['IPTGuM'] = g[3]
        _dbohr_stats['class'] = d['class'].unique()[0]
        fc_stats.append(_dbohr_stats)
fc_stats = pd.concat(fc_stats)

fc_stats.to_csv('../../data/csv/empirical_F_data_statistics.csv', index=False)
