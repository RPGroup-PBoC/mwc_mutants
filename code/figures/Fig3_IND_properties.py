# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.viz
import mut.thermo
import mut.stats
colors = mut.viz.color_selector('mut')
constants = mut.thermo.load_constants()
mut.viz.plotting_style()

# Load the summarized data. 
data = pd.read_csv('../../data/csv/summarized_data.csv')
data = data[(data['class']=='IND')]

# Load the sampling information.
KaKi_samps = pd.read_csv('../../data/csv/KaKi_only_samples.csv')
KaKi_epAI_samps = pd.read_csv('../../data/csv/KaKi_epAI_samples.csv')

# Load the statistics
KaKi_stats = pd.read_csv('../../data/csv/KaKi_only_summary.csv')
KaKi_epAI_stats = pd.read_csv('../../data/csv/DNA_binding_energy_summary.csv')


# Set the figure canvas. 
fig, ax = plt.subplots(2, 3, figsize=(7, 4))

# Set conditional formatting
ax[0, 0].axis('off')
ax[1, 0].axis('off')
ax[0, 1].set_xscale('log')
ax[1, 1].set_xscale('log')

for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=6.5)
    a.yaxis.set_tick_params(labelsize=6.5)
    a.set_ylim([-0.05, 1.2])

# Add appropriate labeling
for i in range(2):
    ax[i, 1].set_ylabel('fold-change', fontsize=6.5)
    ax[i, 1].set_xlabel('IPTG [M]', fontsize=6.5)
    ax[i, 2].set_ylabel('fold-change', fontsize=6.5)
    ax[i, 2].set_xlabel('Bohr parameter [$k_BT$]', fontsize=6.5)

# #############################
# FOLD-CHANGE DATA
# #############################
for g, d in data[data['operator']=='O2'].groupby(['mutant']):
    if g not in ['Q294R', 'Q294K']:
        _ax = ax[0, 1]
    else:
        _ax = ax[1, 1]
            
    _ = _ax.errorbar(d['IPTGuM'] / 1E6, d['mean'], d['sem'],
                    fmt='o', color=colors[g], linestyle='none',
                    linewidth=1, label=g, ms=3)
       
# #############################
# COLLAPSE DATA
# #############################
operator_glyphs = {'O1': 's', 'O2': 'o', 'O3': 'd'}
for g, d in data.groupby(['mutant', 'operator']):
    if g[0] in ['Q294R', 'Q294K']:
        _kaki_samples = KaKi_epAI_samps
        _ax = ax[1, 2]
    else:
        _kaki_samples = KaKi_samps
        _kaki_samples['ep_AI'] = constants['ep_AI']
        _ax = ax[0, 2]
    if (g[1] == 'O3') & (g[0] not in ['Q294K', 'Q294R']):
        op = 'O1'
    else:
        op = g[1]
    # Isolate the mutant and operator
    _params = _kaki_samples[(_kaki_samples['mutant']==g[0]) & (_kaki_samples['operator']==g[1])]
    _mode = _params.iloc[np.argmax(_params['lp__'].values)]
    Ka_samps = _params['Ka']
    Ka_mode = _mode['Ka']
    Ki_samps = _params['Ki']
    Ki_mode = _mode['Ki']
    ep_AI_samps = _params['ep_AI']
    ep_AI_mode  = _mode['ep_AI']
    
    # Instantiate the architecture and compute the mode Bohr parameter. 
    arch = mut.thermo.SimpleRepression(R=d['repressors'].unique(), ep_r=constants[op],
                                      ka=Ka_mode, ki=Ki_mode, ep_ai=ep_AI_mode,
                                      n_sites=constants['n_sites'], n_ns=constants['Nns'],
                                      effector_conc=d['IPTGuM']).bohr_parameter()
    # Plot the collapse data. 
    if g[1] == 'O2':
        face = 'w'
    else:
        face = colors[g[0]]
    _ = _ax.errorbar(arch, d['mean'], d['sem'], linestyle='none',
                    fmt=operator_glyphs[g[1]], color=colors[g[0]],
                    markerfacecolor=face, markeredgecolor=colors[g[0]],
                    markersize=2, linewidth=1, label='__nolegend__')
    

# ###############################
# COLLAPSE CURVE
# ###############################
bohr_range = np.linspace(-8, 8, 200)
master_curve = (1 + np.exp(-bohr_range))**-1
for i in range(2):
    _ = ax[i, 2].plot(bohr_range, master_curve, lw=1, color='k', label='__nolegend__')
    
plt.tight_layout()




