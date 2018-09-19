# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sys.path.insert(0, '../../')
import mut.thermo
import mut.viz
import imp
imp.reload(mut.viz)
pboc = mut.viz.color_selector('pboc')
color = mut.viz.color_selector('mut')
mut.viz.plotting_style()
constants = mut.thermo.load_constants()


# Load the summarized data and trim
data = pd.read_csv('../../data/csv/summarized_data.csv')
data = data[(data['class']=='IND') | (data['class']=='WT')].copy()
data['IPTGM'] = data['IPTGuM'] / 1E6
O2_data = data[data['operator']=='O2']
O1_data = data[data['operator']=='O1']

# Load the necessary samples
kaki_samples = pd.read_csv('../../data/csv/Fig3_O2_KaKi_only_samples.csv')
kaki_stats = pd.read_csv('../../data/csv/Fig3_O2_KaKi_only_stats.csv')
allo_samples = pd.read_csv('../../data/csv/Fig3_O2_KaKi_epAI_samples.csv')
allo_stats = pd.read_csv('../../data/csv/Fig3_O2_KaKi_epAI_stats.csv')

# Define parameter range constants
c_range = np.logspace(-2, 4, 200)
rep_range = np.logspace(0, 4, 200)

# Compute the wild-type curves
wt_induction = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=constants['O2'],
                                          ep_ai=constants['ep_AI'], ka=constants['Ka'],
                                          ki=constants['Ki'], n_ns=constants['Nns'],
                                          effector_conc=c_range).fold_change()
wt_leakiness = mut.thermo.SimpleRepression(R=rep_range, ep_r=constants['O2'],
                                          ep_ai=constants['ep_AI'], ka=constants['Ka'],
                                          ki=constants['Ki'], n_ns=constants['Nns'],
                                          effector_conc=0).fold_change()



# Instantiate the figure.
fig, ax = plt.subplots(2, 2, figsize=(7, 4.5))

# ----------------------------------------------
# SCALING AND LABELS ---------------------------
# ----------------------------------------------
# ax[0,0].axis('off')
ax[0, 1].set_yscale('log')
for a in ax.ravel():
    a.set_xscale('log')
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
  
for a in ax.ravel():
    a.set_ylabel('fold-change', fontsize=8)
ax[0, 1].set_xlabel('repressors per cell', fontsize=8)
ax[0, 0].set_xlabel('IPTG [M]', fontsize=8)
ax[1, 0].set_xlabel('IPTG [M]', fontsize=8)
ax[1, 1].set_xlabel('IPTG [M]', fontsize=8)


# Add titles. 
ax[0,0].set_title('$R$ = 260, O2', fontsize=8, backgroundcolor=pboc['pale_yellow'])
ax[0,1].set_title('$R$ = 260, O2', fontsize=8, backgroundcolor=pboc['pale_yellow'])
ax[1,0].set_title('$R$ = 260, O2', fontsize=8, backgroundcolor=pboc['pale_yellow'])
ax[1,1].set_title('$R$ = 260, O1', fontsize=8, backgroundcolor=pboc['pale_yellow'])
# -----------------------------------------------
# KaKi FITTING ONLY (B) -------------------------
# -----------------------------------------------
for g, d in O2_data.groupby(['mutant']):
    if g == 'wt':
        c='k'
        markerfacecolor = 'w'
        markeredgecolor = 'k'
        ka = constants['Ka']
        ki = constants['Ki']
    else:  
        c = color[g]
        markerfacecolor = c
        markeredgecolor = c 
        ka = kaki_stats[kaki_stats['parameter']=='Ka.{}'.format(g)]['mode'].values[0]
        ki = kaki_stats[kaki_stats['parameter']=='Ki.{}'.format(g)]['mode'].values[0]
    
    
    # Plot the data
    _ = ax[0, 0].errorbar(d['IPTGM'], d['mean'], d['sem'], lw=1, color=c,
                      ms=2, fmt='o', linestyle='none', label=g.upper(), markerfacecolor=markerfacecolor,
                         markeredgecolor=markeredgecolor, markeredgewidth=1)
    
    # Compute the best fit
    if g == 'wt':
        fc = wt_induction

    else:
        fc = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=constants['O2'], ep_ai=constants['ep_AI'],
                                        ka=ka, ki=ki, effector_conc=c_range, n_sites=constants['n_sites'], 
                                        n_ns=constants['Nns']).fold_change()
        cred_region = np.zeros([2, len(c_range)])
        for i, conc in enumerate(c_range):
            arch = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=constants['O2'], ep_ai=constants['ep_AI'],
                                        ka=kaki_samples['Ka.{}'.format(g)], ki=kaki_samples['Ki.{}'.format(g)], 
                                         effector_conc=conc, n_sites=constants['n_sites'], n_ns=constants['Nns']).fold_change()   
            cred_region[:, i] = mut.stats.compute_hpd(arch, 0.95)
            
    # Plot the fits and credible regions
    _ = ax[0, 0].plot(c_range / 1E6, fc, lw=1, color=c, label='__nolegend__')
    
    if g != 'wt':
        _ = ax[0, 0].fill_between(c_range / 1E6, cred_region[0, :], cred_region[1, :], color=c, alpha=0.3, 
                             label='__nolegend__')

# ---------------------------------------------------------------
# LEAKINESS (C) AND CORRECTED INDUCTION PROFILE (D) -------------
# ---------------------------------------------------------------
# Plot the leakiness with WT binding energy and adjusted leakiness
grouped = O2_data[((O2_data['mutant']=='Q294K') |\
                  (O2_data['mutant'] == 'Q294R'))].groupby(['mutant'])

# Plot the theoretical leakiness
_ = ax[0, 1].plot(rep_range, wt_leakiness, 'k-', lw=1, label='WT')

for g, d in grouped:
    # Plot the experimental data for leakiness
    _d = d[d['IPTGuM'] == 0]
    _ = ax[0, 1].errorbar(_d['repressors'], _d['mean'], _d['sem'], color=color[g],
                         lw=1, fmt='o', markersize=3, label=g)
    
    # Plot the experimental data for the induction curve
    _  = ax[1, 0].errorbar(d['IPTGM'], d['mean'], d['sem'], color=color[g],
                       lw=1, fmt='o', markersize=3, label=g)
   
    # Extract the most-likely parameter values
    ep_AI_mode = allo_stats[(allo_stats['parameter']=='ep_AI.{}'.format(g))]['mode'].values[0]
    ka = allo_stats[(allo_stats['parameter']=='Ka.{}'.format(g))]['mode'].values[0]
    ki = allo_stats[(allo_stats['parameter']=='Ki.{}'.format(g))]['mode'].values[0]
    
    # Compute the prediction. 
    fc_leak = mut.thermo.SimpleRepression(R=rep_range, ep_r=constants['O2'],
                                    ka=ka, ki=ki, ep_ai=ep_AI_mode, effector_conc=0, 
                                    n_sites=constants['n_sites'], n_ns=constants['Nns']).fold_change()
    fc_ind = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=constants['O2'],
                                    ka=ka, ki=ki, ep_ai=ep_AI_mode, effector_conc=c_range, 
                                    n_sites=constants['n_sites'], n_ns=constants['Nns']).fold_change()
   # Compute the credible region. 
    leak_cred_region = np.zeros([2, len(rep_range)])
    ind_cred_region = np.zeros([2, len(c_range)])
    for i, r in enumerate(rep_range):
        leak_arch = mut.thermo.SimpleRepression(R=r, ep_r=constants['O2'],
                                        ka=allo_samples['Ka.{}'.format(g)], 
                                        ki=allo_samples['Ki.{}'.format(g)],
                                        ep_ai=allo_samples['ep_AI.{}'.format(g)],
                                        effector_conc=0, n_sites=constants['n_sites'],
                                        n_ns=constants['Nns']).fold_change()
        ind_arch = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=constants['O2'],
                                        ka=allo_samples['Ka.{}'.format(g)], 
                                        ki=allo_samples['Ki.{}'.format(g)],
                                        ep_ai=allo_samples['ep_AI.{}'.format(g)],
                                        effector_conc=c_range[i], n_sites=constants['n_sites'],
                                        n_ns=constants['Nns']).fold_change()

        leak_cred_region[:, i] = mut.stats.compute_hpd(leak_arch, 0.95)
        ind_cred_region[:, i] = mut.stats.compute_hpd(ind_arch, 0.95)
   
    # Plot the predictions and credible regions. 
    _ = ax[0, 1].plot(rep_range, fc_leak, lw=1, color=color[g], label=g)
    _ = ax[0, 1].fill_between(rep_range, leak_cred_region[0, :], leak_cred_region[1, :], color=color[g],
                             alpha=0.3, label='__nolegend__')
    _ = ax[1, 0].plot(c_range / 1E6, fc_ind, lw=1, color=color[g], label=g)
    _ = ax[1, 0].fill_between(c_range / 1E6, ind_cred_region[0, :], ind_cred_region[1, :], color=color[g],
                              alpha=0.3, label='__nolegend__')
    
    
# --------------------------------------------------------------------
# O1 INDUCTION CURVES ------------------------------------------------
# --------------------------------------------------------------------
grouped = O1_data[(O1_data['mutant']=='Q294K') | (O1_data['mutant'] == 'Q294R')].groupby('mutant')
for g, d in grouped:
    # Plot the data. 
    _ = ax[1, 1].errorbar(d['IPTGM'], d['mean'], d['sem'], color=color[g], lw=1, fmt='o', ms=3, label=g)
    
    # Extract the most-likely parameter values 
    ep_AI_mode = allo_stats[(allo_stats['parameter']=='ep_AI.{}'.format(g))]['mode'].values[0]
    ka = allo_stats[(allo_stats['parameter']=='Ka.{}'.format(g))]['mode'].values[0]
    ki = allo_stats[(allo_stats['parameter']=='Ki.{}'.format(g))]['mode'].values[0]
    
    # Compute the best-fit 
    fc = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=constants['O1'],
                                    ka=ka, ki=ki, ep_ai=ep_AI_mode, effector_conc=c_range, 
                                    n_sites=constants['n_sites'], n_ns=constants['Nns']).fold_change()
    
    # Compute the credible regions. 
    cred_region = np.zeros([2, len(c_range)])
    for i, c in enumerate(c_range):
        arch = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=constants['O1'],
                                        ka=allo_samples['Ka.{}'.format(g)], 
                                        ki=allo_samples['Ki.{}'.format(g)],
                                        ep_ai=allo_samples['ep_AI.{}'.format(g)],
                                        effector_conc=c, n_sites=constants['n_sites'],
                                        n_ns=constants['Nns']).fold_change()
        cred_region[:, i] = mut.stats.compute_hpd(arch, 0.95)
        
    # Plot the best fit and credible regions. 
    _ = ax[1, 1].plot(c_range / 1E6, fc, lw=1, color=color[g], label=g)
    _ = ax[1, 1].fill_between(c_range / 1E6, cred_region[0, :], cred_region[1, :], color=color[g], alpha=0.3, label='__nolegend__')

_ = ax[1, 0].legend(fontsize=8, loc='upper left')
_ = ax[1, 1].legend(fontsize=8, loc='upper left')
_ = ax[0, 0].legend(fontsize=8, loc='upper left')
_ = ax[0, 1].legend(fontsize=8, loc='lower left')
plt.tight_layout()
plt.savefig('../figures/Fig3_kaki_plots.pdf', bbox_inches='tight')
