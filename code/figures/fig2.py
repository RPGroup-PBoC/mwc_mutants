# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sys.path.insert(0, '../../')
import mut.viz
import mut.thermo
pboc = mut.viz.color_selector('pboc')
color = mut.viz.color_selector('mut')
induction_colors = sns.color_palette('colorblind', n_colors=7)
induction_colors[4] = sns.xkcd_palette(['dusty purple'])[0]
constants = mut.thermo.load_constants()
mut.viz.plotting_style()

RBS='RBS1027'
# Load the summarized data. 
data = pd.read_csv('../../data/csv/summarized_data.csv')
data = data[(data['operator']=='O2') & 
            ((data['class']=='WT') | (data['class']=='DNA'))].copy()
data['IPTGM'] = data['IPTGuM'] / 1E6
rep_colors = {r:induction_colors[i] for i, r in enumerate(data['repressors'].unique())}

# Load the sampling data. 
epRA_stats = pd.read_csv('../../data/csv/Fig2_O2_DNA_binding_energy_stats.csv')
global_stats = pd.read_csv('../../data/csv/FigS1_O2_DNA_binding_energy_global_stats.csv')

# Compute the wild-type curves. 
rep_range = np.logspace(0, 4, 200)
c_range = np.logspace(-2, 4, 200)
wt_induction = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=constants['O2'],
                                          ep_ai=constants['ep_AI'], ka=constants['Ka'],
                                          ki=constants['Ki'], n_ns=constants['Nns'],
                                          effector_conc=c_range).fold_change()
wt_leakiness = mut.thermo.SimpleRepression(R=rep_range, ep_r=constants['O2'],
                                          ep_ai=constants['ep_AI'], ka=constants['Ka'],
                                          ki=constants['Ki'], n_ns=constants['Nns'],
                                          effector_conc=0).fold_change()

# Instantiate the figure
fig, ax = plt.subplots(3, 2, figsize=(4.5, 6))
_ax = ax.ravel()
mut_ax = {m.upper():_ax[i+2] for i, m in enumerate(data['mutant'].unique())}
               
# Format the axes as necessary
ax[0, 0].axis('off')
ax[0, 1].set_xscale('log')
ax[0, 1].set_yscale('log')
ax[0, 1].set_xlim([1, 1E4])
ax[0, 1].set_ylim([1E-5, 1.2])

# Add sizes and limits
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8) 
for a in ax.ravel()[2:-1]:
    a.set_xscale('log')
    a.set_xlim([1E-8, 1E-2])
    a.set_ylim([-0.05, 1.25])
    a.set_xlabel('IPTG [M]', fontsize=8)
    a.set_ylabel('fold-change', fontsize=8)
    
# Format ticks for parameter estimates
_ax[-1].set_xlim([-0.5, 3.5])
_ax[-1].set_ylim([-16.5, -8])
_ax[-1].set_ylabel('DNA binding energy [$k_BT$]', fontsize=8)
_ax[-1].xaxis.set_ticks([0, 1, 2, 3])
_ax[-1].xaxis.set_ticklabels(np.sort(data['repressors'].unique()).astype(int))

# ------------------------------------
# LEAKINESS (B) ----------------------
# ------------------------------------
# Plot the leakiness
leakiness = data[data['IPTGuM']==0]
for g, d in leakiness.groupby(['mutant']):
    if g.upper() == 'WT':
        # Plot the data
        _ = _ax[1].errorbar(d['repressors'], d['mean'], d['sem'], linestyle='none', fmt='.',
                       ms=3, label='WT', markerfacecolor='w',
                       markeredgecolor='k', markeredgewidth=0.25)
        
        # Plot the theory
        _ = _ax[1].plot(rep_range, wt_leakiness, 'k', lw=1, label='__nolegend__')
    else:
        # Plot the data
        _ = _ax[1].errorbar(d['repressors'], d['mean'], d['sem'], linestyle='none', fmt='.',
                      ms=3, color=color[g], label=g)
        
        # Isolate the parameter values
        epRA = epRA_stats[epRA_stats['parameter']=='ep_RA.{}.{}'.format(g, constants[RBS])].values[0][1:].astype(float)
        
        # Compute the theoretical curves. 
        epRA_mesh, R_mesh = np.meshgrid(epRA, rep_range)
        fc = mut.thermo.SimpleRepression(R=R_mesh, ep_r=epRA_mesh, ep_ai=constants['ep_AI'], effector_conc=0.0,
                                        ka=constants['Ka'], ki=constants['Ki'], n_sites=constants['n_sites'],
                                        n_ns=constants['Nns']).fold_change()
        _ = _ax[1].plot(rep_range, fc.T[0], lw=1, color=color[g])
        
        
# ------------------------------------------------------------- 
# INDUCTION CURVES (C - D)
# -------------------------------------------------------------
for g, d in data.groupby(['mutant', 'repressors']):
    if g[0].upper() == 'WT':
        for a in _ax[2:-1]:
            
            # Plot the data
            _ = a.errorbar(d['IPTGM'], d['mean'], d['sem'], lw=0.75, 
                           linestyle='none', fmt='.', ms=3, 
                           markerfacecolor='w', markeredgecolor='k', 
                           markeredgewidth=0.25, label='WT') 
            # Plot the theory
            _ = a.plot(c_range, wt_induction, 'k', lw=1, label='__nolegend__')
            
    else:
        _ = mut_ax[g[0].upper()].errorbar(d['IPTGM'], d['mean'], 
                                          d['sem'], lw=0.75, linestyle='none', fmt='.',
                                          ms=3, color=rep_colors[g[1]], label=int(g[1]))
        mut_ax[g[0].upper()].set_title(g[0], fontsize=8, backgroundcolor=pboc['pale_yellow'], y=1.04)

# Add the appropriate legend
leg = _ax[2].legend(title='rep. / cell', fontsize=6, handletextpad=0.2, bbox_to_anchor=(0.75, 0.6))
leg.get_title().set_fontsize(6)                                          
                 
# --------------------------------------------------------------
# PARAMETER ESTIMATES (E) 
# --------------------------------------------------------------
# Plot the Points for parameter estimates
rep_idx = {i:r for i, r in enumerate(data['repressors'].unique())}    
for i, r in rep_idx.items():
    for j, dna in enumerate(data['mutant'].unique()):
        # Isolate the parameter estimates. 
        if dna == 'wt':
            _ax[-1].hlines(constants['O2'], -0.5, 3.5, 'k', lw=1, alpha=0.3)
        else:
            param = epRA_stats[epRA_stats['parameter']=='ep_RA.{}.{}'.format(dna, int(r))]
            global_param = global_stats[global_stats['parameter']=='ep_RA.{}'.format(dna)]
            
            _ax[-1].plot(i, param['mode'], '.', ms=3, lw=1, color=color[dna])
            _ax[-1].vlines(i, param['hpd_min'], param['hpd_max'], lw=1, color=color[dna])
            _ax[-1].hlines(global_param['mode'], -0.5, 3.5, color=color[dna], alpha=0.3)
             
plt.tight_layout()
# fin!

