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
model = mut.bayes.StanModel('../stan/empirical_F_data.stan', force_compile=True)


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
    #     _fc_stats = mut.stats.compute_statistics(sample, varnames=['fc_mu'], 
    #                                            logprob_name='lp__')
        _dbohr_stats = mut.stats.compute_statistics(sample, varnames=['delta_bohr'], 
                                               logprob_name='lp__')   
    #     _fc_stats['bohr_mean'] = _bohr_stats['mean']
    #     _fc_stats['bohr_median'] = _bohr_stats['median']
    #     _fc_stats['bohr_mode'] = _bohr_stats['mode']
    #     _fc_stats['bohr_min'] = _bohr_stats['hpd_max']
    #     _fc_stats['bohr_max'] = _bohr_stats['hpd_min']
        _dbohr_stats.rename(columns={'mean':'delta_bohr_mean', 'mode':'delta_bohr_mode',
                                 'median':'delta_bohr_median', 'hpd_min':'delta_bohr_min',
                                 'hpd_max':'delta_bohr_max'}, inplace=True)
        _dbohr_stats['mutant'] = g[0]
        _dbohr_stats['repressors'] = g[1]
        _dbohr_stats['operator'] = g[2]
        _dbohr_stats['IPTGuM'] = g[3]
        _dbohr_stats['class'] = d['class'].unique()[0]
        fc_stats.append(_dbohr_stats)
fc_stats = pd.concat(fc_stats)


# # Sample the posterior
# fit, samples = model.sample(data_dict, iter=5000)

# # Compute the statistics 
# fc_vars = [f'fc_mu[{i}]' for i in data['idx'].unique()]
# bohr_vars = [f'empirical_bohr[{i}]' for i in data['idx'].unique()]
# delta_bohr_vars = [f'delta_bohr[{i}]' for i in data['idx'].unique()]
# fc_stats = mut.stats.compute_statistics(samples, varnames=fc_vars, 
#                                         logprob_name='lp__')
# bohr_stats = mut.stats.compute_statistics(samples, varnames=bohr_vars,
#                                           logprob_name='lp__')
# dbohr_stats = mut.stats.compute_statistics(samples, varnames=bohr_vars,
#                                           logprob_name='lp__')


# # Manually (ugh) merge
# fc_stats['idx'] = data['idx'].unique()
# fc_stats['bohr_mean'] = bohr_stats['mean']
# fc_stats['bohr_median'] = bohr_stats['median']
# fc_stats['bohr_mode'] = bohr_stats['mode']
# fc_stats['bohr_min'] = bohr_stats['hpd_min']
# fc_stats['bohr_max'] = bohr_stats['hpd_max']
# fc_stats['delta_bohr_mean'] = dbohr_stats['mean']
# fc_stats['delta_bohr_median'] = dbohr_stats['median']
# fc_stats['delta_bohr_mode'] = dbohr_stats['mode']
# fc_stats['delta_bohr_min'] = dbohr_stats['hpd_min']
# fc_stats['delta_bohr_max'] = dbohr_stats['hpd_max']

# # Rename the foldchange columns. 
# fc_stats.rename(columns={'mean':'fold_change_mean',
#                         'median':'fold_change_median',
#                         'hpd_min':'fold_change_min',
#                         'hpd_max': 'fold_change_max'},
#                inplace=True)


# # Add identifiers
# for i in data['idx'].unique():
#     fc_stats.loc[fc_stats['idx']==i,
#                  'mutant'] = data[data['idx']==i]['mutant'].unique()
#     fc_stats.loc[fc_stats['idx']==i,
#                  'IPTGuM'] = data[data['idx']==i]['IPTGuM'].unique()
#     fc_stats.loc[fc_stats['idx']==i,
#                  'repressors'] = data[data['idx']==i]['repressors'].unique()
#     fc_stats.loc[fc_stats['idx']==i,
#                  'operator'] = data[data['idx']==i]['operator'].unique()
#     fc_stats.loc[fc_stats['idx']==i,
#                  'class'] = data[data['idx']==i]['class'].unique()

fc_stats.to_csv('../../data/csv/empirical_F_statistics.csv', index=False)

# Compute the posterior predictive checks
# rep_vars = [f'y_rep[{i+1}]' for i in range(len(data))]
# rep_stats = mut.stats.compute_statistics(samples, varnames=rep_vars, 
#                                          logprob_name='lp__')
# rep_stats.to_csv('../../data/csv/inferred_fc_ppc.csv', index=False)


