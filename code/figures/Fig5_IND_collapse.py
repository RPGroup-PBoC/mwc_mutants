# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import mut.thermo
import mut.stats
import mut.viz
import seaborn as sns
_colors = sns.color_palette('deep')
op_colors = {'O1':_colors[0], 'O2':_colors[1], 'O3':_colors[2]}
pboc = mut.viz.color_selector('pboc')
colors = mut.viz.color_selector('mut')
constants = mut.thermo.load_constants()
mut.viz.plotting_style()

# Load the data
data = pd.read_csv('../../data/csv/compiled_data.csv')
ind_data = data[data['class']=='IND'].copy()
kaki_only_stats = pd.read_csv('../../data/csv/KaKi_only_summary.csv')
kaki_epAI_stats = pd.read_csv('../../data/csv/KaKi_epAI_summary.csv')
kaki_epAI_samps = pd.read_csv('../../data/csv/KaKi_epAI_samples.csv')
delta_F = pd.read_csv('../../data/csv/empirical_F_statistics.csv')
ind_delF = delta_F[(delta_F['class'] == 'IND')]

# Go through and compute the bohr parameter
for g, d in ind_data.groupby(['operator', 'mutant']):
    # Get the statistics. 
    _stats = kaki_epAI_stats[(kaki_epAI_stats['mutant']==g[1]) & 
                             (kaki_epAI_stats['operator']=='O2')]
    _samps = kaki_epAI_samps[(kaki_epAI_samps['mutant']==g[1]) & 
                        (kaki_epAI_samps['operator']=='O2')]
    ka_median = _stats[_stats['parameter']=='Ka']['median'].values[0]
    ki_median = _stats[_stats['parameter']=='Ki']['median'].values[0]
    epAI_median =  _stats[_stats['parameter']=='ep_AI']['median'].values[0]
    median_bohr = mut.thermo.SimpleRepression(R=260, ep_r=constants[g[0]], 
            ka=ka_median, ki=ki_median, ep_ai=epAI_median, 
            n_sites=constants['n_sites'], effector_conc=d['IPTGuM']).bohr_parameter()
    ind_data.loc[(ind_data['mutant']==g[1]) & 
    (ind_data['operator']==g[0]), 'bohr'] = median_bohr
    cred_region = np.zeros((2, len(d)))
    for i, c in enumerate(d['IPTGuM'].values):
        arch = mut.thermo.SimpleRepression(R=260, ep_r=constants[g[0]], 
                ka=_samps['Ka'], ki=_samps['Ki'], n_sites=constants['n_sites'],
                effector_conc=c, ep_ai=_samps['ep_AI']).bohr_parameter()
        cred_region[:, i] = mut.stats.compute_hpd(arch, 0.95)
    ind_data.loc[(ind_data['mutant']==g[1]) & 
    (ind_data['operator']==g[0]), 'bohr_min'] = cred_region[0, :]
    ind_data.loc[(ind_data['mutant']==g[1]) & 
    (ind_data['operator']==g[0]), 'bohr_max'] = cred_region[1,: ]
ind_summary = ind_data.groupby(['mutant', 'operator', 'IPTGuM']).agg(('mean', 'sem')).reset_index()

# Define the  constants for plotting the fits. 
c_range = np.logspace(-3, 4, 500)
c_range[0] = 0
bohr_range = np.linspace(-8, 8, 500)

# Set up the figure canvas
fig, ax = plt.subplots(2, 4, figsize=(7, 4), sharey=True)

# Add labels and format the axes
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_ylim([-0.1, 1.2])
for i in range(4):
    ax[0, i].set_xscale('symlog', linthreshx=1E-3, linscalex=0.5)
    ax[0, i].set_xlabel('IPTG [µM]', fontsize=8)
    ax[1, i].set_xlabel('free energy [$k_BT$]', fontsize=8)
    ax[1, i].set_xlim([-8, 8])
    ax[0, i].set_xlim([0, 1E4])
    ax[0, i].set_xlim([-0.001, 1E4])

mut_ind = {'F164T':3, 'Q294V':2, 'Q294K':1, 'Q294R':0}
for m, i in mut_ind.items():
    ax[0,i].set_title(m, fontsize=8, backgroundcolor=pboc['pale_yellow'], y=1.08)
    ax[0, i].set_xticks([0, 1E-2, 1E0, 1E2, 1E4])

ax[0, 0].set_ylabel('fold-change', fontsize=8)
ax[1, 0].set_ylabel('fold-change', fontsize=8)
# ##############################################################################
#  INDUCTION AND COLLAPSE DATA
# ##############################################################################
for g, d in ind_summary.groupby(['mutant', 'operator']):    
    fc_ax = ax[0, mut_ind[g[0]]]
    bohr_ax = ax[1, mut_ind[g[0]]]
    if g[1] == 'O2':
        face = 'w'
        zorder = 1000
    else:
        face = op_colors[g[1]]
        zorder = 1

    label = g[1]
    fc_ax.errorbar(d['IPTGuM'], d['fold_change']['mean'], 
                   d['fold_change']['sem'], fmt='.', ms=5, linestyle='none',
                    lw=1, capsize=1,  markerfacecolor=face, 
                    color=op_colors[g[1]], zorder=zorder, label=label)
    bohr_ax.errorbar(d['bohr']['mean'], d['fold_change']['mean'], 
                   d['fold_change']['sem'], fmt='.', ms=5, linestyle='none',
                    lw=1, capsize=1,
                    markerfacecolor=face, color=op_colors[g[1]],
                    zorder=zorder, label=label)

# ##############################################################################
# INDUCTION CURVES  
# ##############################################################################
for g, d in ind_summary.groupby(['mutant', 'operator']):
    _axis = ax[0, mut_ind[g[0]]]

    # Plot the "poor" prediction curves
    _stats = kaki_only_stats[(kaki_only_stats['mutant']==g[0]) & 
                            (kaki_only_stats['operator']=='O2')]
    ka_median = _stats[_stats['parameter']=='Ka']['median'].values[0]
    ki_median = _stats[_stats['parameter']=='Ki']['median'].values[0]
    fit = mut.thermo.SimpleRepression(R=260, ep_r=constants[g[1]],
                                    ka=ka_median, ki=ki_median, 
                                    ep_ai=constants['ep_AI'], 
                                    n_sites=constants['n_sites'], 
                                    effector_conc=c_range).fold_change()
    _axis.plot(c_range, fit, color=op_colors[g[1]], lw=1, label='__nolegend__',
            linestyle=':')


    # Plot the "good" fits
    _stats = kaki_epAI_stats[(kaki_epAI_stats['mutant']==g[0]) & 
                            (kaki_epAI_stats['operator']=='O2')]
    ka_median = _stats[_stats['parameter']=='Ka']['median'].values[0]
    ki_median = _stats[_stats['parameter']=='Ki']['median'].values[0]
    epAI_median = _stats[_stats['parameter']=='ep_AI']['median'].values[0]
    _samps = kaki_epAI_samps[(kaki_epAI_samps['mutant']==g[0]) &
                             (kaki_epAI_samps['operator']=='O2')]
    # fit = mut.thermo.SimpleRepression(R=260, ep_r=constants[g[1]],
                                    # ka=ka_median, ki=ki_median, 
                                    # ep_ai=epAI_median, 
                                    # n_sites=constants['n_sites'], 
                                    # effector_conc=c_range).fold_change()
    # _axis.plot(c_range, fit, color=op_colors[g[1]], lw=1, label='__nolegend__',
            # linestyle=':')
    cred_region = np.zeros((2, len(c_range)))
    for i, c in enumerate(c_range):
        arch = mut.thermo.SimpleRepression(R=260, ep_r=constants[g[1]],
                                    ka=_samps['Ka'], ki=_samps['Ki'], 
                                    ep_ai=_samps['ep_AI'], 
                                    n_sites=constants['n_sites'], 
                                    effector_conc=c).fold_change()
        cred_region[:, i] = mut.stats.compute_hpd(arch, 0.95)

    

    _axis.fill_between(c_range, cred_region[0, :], cred_region[1, :], color=op_colors[g[1]],
    alpha=0.5)

# ##############################################################################
# COLLAPSE MASTER CURVE
# ##############################################################################
for i in range(4):
    ax[1, i].plot(bohr_range, (1 + np.exp(-bohr_range))**-1, 'k-', lw=1,
                    zorder=1)

# ##############################################################################
# LEGEND
# ##############################################################################
leg = ax[0, 0].legend(title=r'operator', fontsize=7,
    handletextpad=0.1)
leg.get_title().set_fontsize(7)
plt.tight_layout()
plt.savefig('../../figures/Fig5_ind_collapse.pdf', bbox_inches='tight')
