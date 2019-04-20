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
model = mut.bayes.StanModel('../stan/empirical_F.stan') #, force_compile=True)

# Define the parameters
n_rep = 10
n_points = 200

# Generate the fake datasets
F_mu = np.linspace(-8, 8, n_points)

fc_mu2 = (1 + np.exp(-F_mu))**-1
fc_mu = (1 + np.exp(-F_mu))**-1

# Draw sigmas out of a half normal 
sig = np.abs(np.random.normal(0, 0.1, len(F_mu)))
dfs = []
for i in range(n_points):
    fc_rand = np.random.normal(fc_mu[i], sig[i], n_rep)
    df = pd.DataFrame([])
    df['fold_change'] = fc_rand
    df['fc_mu'] = fc_mu[i]
    df['fc_mu_ref'] = fc_mu2[i]
    df['fc_sig'] = sig[i]
    df['bohr']  = F_mu[i]
    df['draw'] = i
    dfs.append(df)
df = pd.concat(dfs)

samp_dfs = []
stat_dfs = []
for g, d in tqdm.tqdm(df.groupby(['draw'])):

    # Sample the posterior for the data
    _, samples = model.sample(dict(N=len(d), foldchange=d['fold_change']), 
                            control=dict(adapt_delta=0.95))
    samples['true_mu'] = d['fc_mu'].values[0]
    samples['true_sig'] = d['fc_sig'].values[0]
    samples['true_bohr'] = d['bohr'].values[0]
    samp_dfs.append(samples)

    # Sample fake data around the mean with the average sigma to estimate error
    sim_data = np.random.normal(d['fc_mu_ref'].values[0], samples['fc_sigma'].median(), 
                              len(d))
    _, samples_sim = model.sample(dict(N=len(sim_data),foldchange=sim_data),
                         control=dict(adapt_delta=0.95))

    # Compute the empirical bohr and delta F 
    samples['sim_empirical_bohr'] = -np.log(samples_sim['fc_mu']**-1 -1)
    samples['empirical_bohr'] = -np.log(samples['fc_mu']**-1 - 1)
    samples['correction'] = d['bohr'].values[0] - samples['sim_empirical_bohr'].median()
    samples['delta_bohr'] = d['bohr'].values[0] - samples['empirical_bohr'] 
    samples['ref_delta_bohr'] = d['bohr'].values[0] - samples['sim_empirical_bohr']
    samples['corr_delta_bohr'] = samples['sim_empirical_bohr']  - samples['empirical_bohr']



    # Identify the extrema
    extrema = (samples['fc_mu'] < samples['fc_sigma']).astype(int) +\
              (1-samples['fc_mu'] < samples['fc_sigma']).astype(int)

    # Compute the delta F error of the reference, given the sigma
    delta_F_ref_upper = np.nan_to_num(d['bohr'].values[0] +\
         np.log((samples['true_mu'] + samples['fc_sigma'])**-1 - 1))
    delta_F_ref_lower = np.nan_to_num(d['bohr'].values[0] +\
         np.log((samples['true_mu'] -samples['fc_sigma'])**-1 - 1))
    samples['correction'] = (delta_F_ref_upper - delta_F_ref_lower) * extrema 
    samples['delta_bohr_corrected'] = samples['delta_bohr'] +\
                                            np.sign(d['bohr'].values[0]) * samples['correction']

    _dbohr_stats = mut.stats.compute_statistics(samples, varnames=['delta_bohr', 'empirical_bohr', 'fc_mu', 
                                                'fc_sigma', 'delta_bohr_corrected', 
                                                'correction', 
                                                'sim_empirical_bohr', 
                                                'corr_delta_bohr', 
                                                'ref_delta_bohr'], 
                                               logprob_name='lp__')   
    _dbohr_stats['true_mu'] = d['fc_mu'].values[0] 
    _dbohr_stats['true_sig'] = d['fc_sig'].values[0]
    _dbohr_stats['true_bohr'] = d['bohr'].values[0]
    stat_dfs.append(_dbohr_stats)

samples = pd.concat(samp_dfs)
stats = pd.concat(stat_dfs)
df.to_csv('../../data/csv/pathological_F_data.csv', index=False)
samples.to_csv('../../data/csv/pathological_F_samples.csv', index=False)
stats.to_csv('../../data/csv/pathological_F_stats.csv', index=False)



plt.semilogx(fc_mu, F_mu) 