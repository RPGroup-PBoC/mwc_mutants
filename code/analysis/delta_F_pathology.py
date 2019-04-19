# -*- coding: utf-8 -*-
import sys
sys.path.insert(0,'../../')
import numpy as np
import pandas as pd
import mut.thermo
import mut.bayes
import tqdm
constants = mut.thermo.load_constants()

# Load the statistical model. 
model = mut.bayes.StanModel('../stan/empirical_F.stan')

# Define the parameters
n_rep = 10
n_points = 200
# Generate the fake dataset
c_range = np.logspace(-4, 6, n_points)
ep_r = -16
ka = 200
ep_AI = 5
ki = 0.1
R = 100
arch = mut.thermo.SimpleRepression(R, ep_r, ka=ka, ki=ki, ep_ai=ep_AI,
                                   effector_conc=c_range)
fc_mu = arch.fold_change()
F_mu = arch.bohr_parameter()

# Draw sigmas out of a half normal 
sig = np.abs(np.random.normal(0, 0.05, len(c_range)))
dfs = []
for i in range(n_points):
    fc_rand = np.random.normal(fc_mu[i], sig[i], n_rep)
    df = pd.DataFrame([])
    df['fold_change'] = fc_rand
    df['fc_mu'] = fc_mu[i]
    df['fc_sig'] = sig[i]
    df['bohr']  = F_mu[i]
    df['draw'] = i
    df['IPTGuM'] = c_range[i]
    dfs.append(df)
df = pd.concat(dfs)

samp_dfs = []
stat_dfs = []
for g, d in tqdm.tqdm(df.groupby(['IPTGuM'])):
    _, samples = model.sample(dict(N=len(d), foldchange=d['fold_change']), 
                            control=dict(adapt_delta=0.95))
    samples['IPTGuM'] = g
    samples['true_mu'] = d['fc_mu'].values[0]
    samples['true_sig'] = d['fc_sig'].values[0]
    samples['true_bohr'] = d['bohr'].values[0]
    samp_dfs.append(samples)

    # Compute the empirical bohr and delta F 
    samples['empirical_bohr'] = -np.log(samples['fc_mu']**-1 - 1)
    samples['delta_bohr'] = d['bohr'].values[0] - samples['empirical_bohr']

    # Identify the extrema
    extrema = (samples['fc_mu'] < samples['fc_sigma']).astype(int) +\
              (1-samples['fc_mu'] < samples['fc_sigma']).astype(int)

    # Compute the delta F error of the reference, given the sigma
    delta_F_ref_upper = np.nan_to_num(d['bohr'].values[0] +\
         np.log((samples['true_mu'].values + samples['fc_sigma'])**-1 - 1))
    delta_F_ref_lower = np.nan_to_num(d['bohr'].values[0] +\
         np.log((samples['true_mu'].values - samples['fc_sigma'])**-1 - 1))
    samples['correction'] = (delta_F_ref_upper-delta_F_ref_lower) * extrema
    samples['delta_bohr_corrected'] = samples['delta_bohr'] +\
                                            np.sign(d['bohr'].values[0]) * samples['correction']

    _dbohr_stats = mut.stats.compute_statistics(samples, varnames=['delta_bohr', 'empirical_bohr', 'fc_mu', 
                                                'fc_sigma', 'delta_bohr_corrected', 'correction'], 
                                               logprob_name='lp__')   
    _dbohr_stats['IPTGuM'] = g
    _dbohr_stats['true_mu'] = d['fc_mu'].values[0] 
    _dbohr_stats['true_sig'] = d['fc_sig'].values[0]
    _dbohr_stats['true_bohr'] = d['bohr'].values[0]
    stat_dfs.append(_dbohr_stats)

samples = pd.concat(samp_dfs)
stats = pd.concat(stat_dfs)
df.to_csv('../../data/csv/pathological_F_data.csv', index=False)
samples.to_csv('../../data/csv/pathological_F_samples.csv', index=False)
stats.to_csv('../../data/csv/pathological_F_stats.csv', index=False)