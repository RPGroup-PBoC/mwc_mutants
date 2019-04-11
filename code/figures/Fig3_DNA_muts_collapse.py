# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
import mut.stats
import seaborn as sns
pboc = mut.viz.color_selector('pboc')
colors = mut.viz.color_selector('mut')
constants = mut.thermo.load_constants()
mut.viz.plotting_style()
fit_strain = 260 # In repressors per cell
rep_colors = {60:pboc['blue'], 124:pboc['purple'], 260:pboc['red'], 1220:pboc['green']}

# ########################################################
# DATA LOADING & PRUNING
# ########################################################
fc_data = pd.read_csv('../../data/csv/compiled_data.csv') 
epRA_stats = pd.read_csv('../../data/csv/DNA_binding_energy_summary.csv')
DNA_fc = fc_data[((fc_data['class']=='DNA') | (fc_data['class']=='WT')) & 
                (fc_data['operator']=='O2')]
grouped = DNA_fc.groupby(['mutant', 'repressors', 
                          'IPTGuM']).agg(('mean', 'sem')).reset_index()
c_range = np.logspace(-2, 4, 500)
c_range[0] = 0 
bohr_range = np.linspace(-8, 8, 500)


# ########################################################
# FIGURE INSTANTIATION / FORMATTING 
# ########################################################
fig, ax = plt.subplots(2, 3, figsize=(5, 3.5), sharey=True)
axes = {'Y20I': 0, 'Q21A':1, 'Q21M':2}

# Format axes
for i in range(3): 
    ax[0, i].set_ylim([-0.1, 1.2])
    ax[1, i].set_xlim([-8, 8])
    ax[0, i].set_xscale('symlog', linthreshx=1E-2)
    ax[0, i].set_xticks([0, 1E-2, 1E0, 1E2, 1E4])
    ax[0, i].set_xlabel('IPTG  [µM]', fontsize=8)
    ax[0, i].set_xlim([-0.001, 1E4])
    ax[1, i].set_xlabel('free energy [$k_BT$]', fontsize=8)
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8) 

# Add labels
for m, a in axes.items():
    ax[0, a].set_title(m, y=1.08, backgroundcolor=pboc['pale_yellow'],
                      fontsize=8)
for i in range(2):
    ax[i, 0].set_ylabel('fold-change', fontsize=8)

# ######################################################################
# COLLAPSE CURVE
# ######################################################################
for i in range(3):
    _ = ax[1, i].plot(bohr_range, (1 + np.exp(-bohr_range))**-1, lw=0.75, color='k', linestyle='-')
    
 
# #######################################
# INDUCTION DATA 
# ######################################
for g, d in grouped[grouped['mutant'] != 'wt'].groupby(['mutant', 'repressors']): 
    if g[1] == fit_strain:
        face = 'w'
        zorder=1000
    else:
        face = rep_colors[g[1]]
        zorder=1
    _ax = ax[0, axes[g[0]]]
    _ax.errorbar(d['IPTGuM'], d['fold_change']['mean'], d['fold_change']['sem'], lw=0.75, 
                capsize=1, linestyle='none', fmt='.', ms=5, markerfacecolor=face,
                color=rep_colors[g[1]], zorder=zorder, label=int(g[1]))
    
    
leg = ax[0, -1].legend(loc='upper left', fontsize=6, title='rep / cell')
leg.get_title().set_fontsize(6)

    
# #############################################################################
# MUTANT PROFILE CURVES 
# #############################################################################
fit_stats = epRA_stats[(epRA_stats['repressors']==fit_strain) & (epRA_stats['parameter']=='ep_RA')]

for g, d in fit_stats[fit_stats['mutant'] != 'wt'].groupby(['mutant']):
    _ind_ax = ax[0, axes[g]]
    for r, c in rep_colors.items():
        _c, _ep = np.meshgrid(c_range, d[['median', 'hpd_min', 'hpd_max']].values)
        arch = mut.thermo.SimpleRepression(R=r, ep_r=_ep, 
                                           ka=constants['Ka'],
                                         ki=constants['Ki'], ep_ai=constants['ep_AI'],
                                         n_sites=constants['n_sites'], n_ns=constants['Nns'],
                                         effector_conc=_c)
        _fc = arch.fold_change()
        _ind_ax.fill_between(c_range, _fc[1, :], _fc[2, :], color=c, alpha=0.5)

        
# ###############################################################################
# MUTANT COLLAPSE CURVES
# ################################################################################
for g, d in fit_stats[fit_stats['mutant'] != 'wt'].groupby('mutant'):
    _bohr_ax = ax[1, axes[g]]
    for r, c in rep_colors.items():
        _mut_strain = grouped[(grouped['mutant']==g) & (grouped['repressors']==r)]
        if r == 260:
            face='w'
            zorder=1000
        else:
            face = c
            zorder=100
        _c, _ep = np.meshgrid(_mut_strain['IPTGuM'], d[['median', 'hpd_min', 'hpd_max']])
        arch = mut.thermo.SimpleRepression(R=r, ep_r=_ep, ka=constants['Ka'],
                                          ki=constants['Ki'], ep_ai=constants['ep_AI'],
                                          n_sites=constants['n_sites'], effector_conc=_c)
        bohr = arch.bohr_parameter()
        _bohr_ax.errorbar(bohr[0, :], _mut_strain['fold_change']['mean'], 
                          _mut_strain['fold_change']['sem'], ms=5, markerfacecolor=face,
                          color=c, linestyle='none', fmt='.', zorder=zorder)
        
        _bohr_ax.hlines(_mut_strain['fold_change']['mean'], bohr[1, :], bohr[2, :], color=c,
                       linewidth=1, zorder=zorder)
    
plt.tight_layout()
plt.savefig('../../figures/Fig3_DNA_collapse.pdf', bbox_inches='tight')
