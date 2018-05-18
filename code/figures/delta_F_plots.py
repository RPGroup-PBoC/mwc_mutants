# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, '../../')
import mut.thermo
import mut.stats
import mut.viz 
mut.viz.plotting_style()
pboc_colors = mut.viz.color_selector('pboc')

data = pd.read_csv('../../data/csv/compiled_data.csv')
muts = data[(data['mutant'] != 'wt') & (data[data['mutant'] != 'Q294R'])]['mutant'].unique()

# Load the chains.
DNA_chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_DNA_strict.csv')
IND_chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_IND_strict.csv')
DBL_chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_DBL_strict.csv')
global_chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_global_strict.csv')

#%% Compute the statistics
DNA_stats = mut.stats.compute_statistics(DNA_chains, logprob_name='lnprobability')
IND_stats = mut.stats.compute_statistics(IND_chains, logprob_name='lnprobability')
DBL_stats = mut.stats.compute_statistics(DBL_chains, logprob_name='lnprobability')
global_stats = mut.stats.compute_statistics(global_chains, logprob_name='lnprobability')
stats = pd.concat([DNA_stats, IND_stats], ignore_index=True)
global_stats = pd.concat([global_stats, DBL_stats], ignore_index=True)

# %%
# Compute the WT bohr parameter.
wt_eps_r = -13.9
wt_ka = 139E-6
wt_ki = 0.53E-6
wt_bohr = mut.thermo.SimpleRepression(R=260, ep_r=wt_eps_r, ka=wt_ka, ki=wt_ki, ep_ai=4.5,
                                       effector_conc=1E9).bohr_parameter()
pred_delta = {}
pred_err = {}
meas_delta = {}
meas_err = {}
# Loop through each mutant and calculated the predicted F.
for i, m in enumerate(muts):
    if m in ['Y20I', 'Q21A', 'Q21M']:
        # Compute the predicted bohr parameter.
        pred_ep_r = stats[stats['parameter']=='{}_eps_r'.format(m)]['mode'].values[0]
        pred_ka = wt_ka
        pred_ki = wt_ki

        # isolate the chain for an HPD
        chains_ep_r = DNA_chains['{}_eps_r'.format(m)]
        chains_ka = wt_ka
        chains_ki = wt_ki
    elif m in ['F164T', 'Q294V']:
        pred_ep_r = wt_eps_r
        pred_ka = np.exp(stats[stats['parameter']=='{}_ka'.format(m)]['mode'].values[0])
        pred_ki = np.exp(stats[stats['parameter']=='{}_ka'.format(m)]['mode'].values[0])
        chains_ep_r = wt_eps_r
        chains_ka = np.exp(IND_chains['{}_ka'.format(m)])
        chains_ki = np.exp(IND_chains['{}_ki'.format(m)])
        
    
    elif 'Q294K' not in m:
        DNA_mut = m.split('-')[0]
        IND_mut = m.split('-')[1]
        if IND_mut != 'Q294K':
            pred_ep_r = stats[stats['parameter']=='{}_eps_r'.format(DNA_mut)]['mode'].values[0]
            pred_ka = np.exp(stats[stats['parameter']=='{}_ka'.format(IND_mut)]['mode'].values[0])
            pred_ki = np.exp(stats[stats['parameter']=='{}_ka'.format(IND_mut)]['mode'].values[0])
            # chains_ep_r = DNA_chains['{}_eps_r'.format(DNA_mut)]
            # chains_ka = np.exp(IND_chains['{}_ka'.format(IND_mut)])
            # chains_ki = np.exp(IND_chains['{}_ki'.format(IND_mut)])
 
    # Compute the measured parameter
    meas_ep_r = global_stats[global_stats['parameter']=='{}_eps'.format(m)]['mode'].values[0]
    meas_ka = np.exp(global_stats[global_stats['parameter']=='{}_ka'.format(m)]['mode'].values[0])
    meas_ki = np.exp(global_stats[global_stats['parameter']=='{}_ki'.format(m)]['mode'].values[0])

    if len(m.split('-')) == 2:
        meas_chains_ep_r = DBL_chains['{}_eps'.format(m)]
        meas_chains_ka = np.exp(DBL_chains['{}_ka'.format(m)])
        meas_chains_ki = np.exp(DBL_chains['{}_ki'.format(m)])
    else:
        meas_chains_ep_r = global_chains['{}_eps'.format(m)]
        meas_chains_ka = np.exp(global_chains['{}_ka'.format(m)])
        meas_chains_ki = np.exp(global_chains['{}_ki'.format(m)])


    # Assemble the architectures and compute.  
    predicted_bohr = mut.thermo.SimpleRepression(R=260, ep_r=pred_ep_r, ka=pred_ka/1E6, ki=pred_ki/1E6, ep_ai=4.5,
                                              effector_conc=1E9).bohr_parameter()

    measured_bohr = mut.thermo.SimpleRepression(R=260, ep_r=meas_ep_r, ka=meas_ka/1E6, ki=meas_ki/1E6, ep_ai=4.5,
                                              effector_conc=1E9).bohr_parameter()
    measured_bohr_err = mut.stats.compute_hpd(mut.thermo.SimpleRepression(R=260, ep_r=meas_chains_ep_r, ka=meas_chains_ka/1E6, ki=meas_chains_ki/1E6, ep_ai=4.5,
                                              effector_conc=1E9).bohr_parameter(), 0.95)
    pred_err[m] = (wt_bohr - predicted_bohr_err[0], wt_bohr - predicted_bohr_err[1])
    meas_err[m] = (wt_bohr - measured_bohr_err[0], wt_bohr - measured_bohr_err[1])
    pred_delta[m] = wt_bohr - predicted_bohr
    meas_delta[m] = wt_bohr - measured_bohr

# %%
# Define the plotting order
order = ['Y20I', 'Q21A', 'Q21M', 'F164T', 'Q294V', 'Y20I-F164T', 
        'Y20I-Q294V', 'Q21A-F164T', 'Q21A-Q294V', 'Q21M-F164T', 
        'Q21M-Q294V']
index = {o:i for i, o in enumerate(order)}

fig, ax = plt.subplots(1, 1, figsize=(4, 6))
ax.tick_params(labelsize=8)

for delta in pred_delta:
    if ('Q294K' not in delta) & ('Q294R' not in delta):
        _ = ax.plot(pred_delta[delta], index[delta] + 0.2, 'o', color=pboc_colors['red'])
        _ = ax.hlines(index[delta] + 0.2, pred_err[m], color=pboc_colors['red'])
        _ = ax.plot(meas_delta[delta], index[delta] - 0.2, 'o', color=pboc_colors['blue'])
        _ = ax.hlines(index[delta] - 0.2, meas_err[m], color=pboc_colors['red'])

ax.set_xlabel('$\Delta F$  $(k_BT)$', fontsize=8)
ax.set_ylabel('mutant', fontsize=8)
ax.set_yticks(np.arange(0, len(order), 1))
ax.set_yticklabels(order)

# %%
m