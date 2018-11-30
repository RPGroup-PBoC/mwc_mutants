# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mut.stats
import mut.thermo
import mut.viz
constants = mut.thermo.load_constants()
colors = mut.viz.color_selector('mut')
mut.viz.plotting_style()

# Load the data
data = pd.read_csv('../../data/csv/summarized_data.csv')

# Restrict the data set to Q294R and Q294K
data = data[(data['mutant']=='Q294R') | (data['mutant']=='Q294K')]

# Load the samples. 
KaKi_R_samples = pd.read_csv('../../data/csv/KaKi_R_samples.csv')
KaKi_R_stats = pd.read_csv('../../data/csv/KaKi_R_summary.csv')

# Define the ranges for plotting of values
c_range = np.logspace(-4, 2, 200)
bohr_range = np.linspace(-10, 10, 200)

# Instantiate the figure
fig = plt.figure(figsize=(7, 3.5))
gs = gridspec.GridSpec(3, 6)

# Set up the axes
ax = []
for i in range(3):
    _ax_row = []
    for j in range(3):
        _ax = fig.add_subplot(gs[i, j])
        _ax_row.append(_ax)
    ax.append(_ax_row)
ax = np.array(ax)    
collapse_ax = fig.add_subplot(gs[:, 3:])     

# Add appropriate scaling
for a in ax.ravel():
    a.set_xscale('log')
    a.set_ylim([-0.05, 1.2])
    a.set_xlim([1E-8, 1E-2])
    
# Format the axes. 
for i in range(3):
    for j in range(3):
        if i == j:
            ax[i, j].set_facecolor('#ecf0f1')
        if i < 2:
            ax[i, j].xaxis.set_ticklabels([])
        if j > 0:
            ax[i, j].yaxis.set_ticklabels([])
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=6.5)
    a.yaxis.set_tick_params(labelsize=6.5) 
collapse_ax.xaxis.set_tick_params(labelsize=6.5)
collapse_ax.yaxis.set_tick_params(labelsize=6.5)

# ##############################
# INDUCTION DATA 
# ##############################
# Set the indices for the operators. 
op_ax = {'O1':0, 'O2':1, 'O3':2}

for g, d in data.groupby(['mutant', 'operator']):
    # Identify the axes. 
    for i in range(3):
        if i == op_ax[g[1]]:
            face = 'w'
        else:
            face = colors[g[0]]
        ax[op_ax[g[1]], i].errorbar(d['IPTGuM'] / 1E6, d['mean'], d['sem'],
                                fmt='o', color=colors[g[0]], lw=1,
                                markerfacecolor=face, ms=3, label=g[0])
        
# #############################
# INDUCTION PROFILES
# #############################
for i, m in enumerate(['Q294R', 'Q294K']):
    for j, op_pred in enumerate(op_ax.keys()):
        # Isolate the statistics. 
        stats = KaKi_R_stats[(KaKi_R_stats['mutant']==m) & 
                            (KaKi_R_stats['operator']==op_pred)]
        samps = KaKi_R_samples[(KaKi_R_samples['mutant']==m) &
                              (KaKi_R_samples['operator']==op_pred)]
        Ka = stats[stats['parameter']=='Ka']['mode'].values[0]
        Ka_samps = samps['Ka']
        Ki = stats[stats['parameter']=='Ki']['mode'].values[0]
        Ki_samps = samps['Ki']
        R = stats[stats['parameter']=='R']['mode'].values[0]
        R_samps = samps['R']
        for k, op_fit in enumerate(op_ax.keys()):  
            # Compute the best-fit curve
            arch = mut.thermo.SimpleRepression(R=R, ep_r=constants[op_fit], effector_conc=d['IPTGuM'],
                                       ka=Ka, ki=Ki, ep_ai=constants['ep_AI'], n_ns=constants['Nns'],
                                      n_sites=constants['n_sites']).fold_change()
            
            # Compute the credible region
            cred_region = np.zeros((2, len(c_range)))
            for z, c in enumerate(c_range):
                _arch = mut.thermo.SimpleRepression(R=R_samps, ep_r=constants[op_fit], effector_conc=d['IPTGuM'],
                                       ka=Ka_samps, ki=Ki_samps, ep_ai=constants['ep_AI'], n_ns=constants['Nns'],
                                      n_sites=constants['n_sites']).fold_change()
                cred_region[:, z] = mut.stats.compute_hpd(_arch, 0.95)
                
            # Plot the line of best fit and the credible regions
            _ = ax[j, k].plot(c_range / 1E6, arch, color=colors[m], lw=1, legend='__nolegend__')
            _ = ax[j, k].fill_between(c_range / 1E6, cred_region[0, :], cred_region[1, :],
                                      color=colors[m], lw=1, legend='__nolegend__')
             
# #############################
# COLLAPSE DATA 
# #############################
for g, d in data.groupby(['mutant', 'operator']):
    # Isolate the statistics. 
    stats = KaKi_R_stats[(KaKi_R_stats['mutant']==g[0]) & 
                        (KaKi_R_stats['operator']=='O2')]
    Ka = stats[stats['parameter']=='Ka']['mode'].values[0]
    Ki = stats[stats['parameter']=='Ki']['mode'].values[0]
    R = stats[stats['parameter']=='R']['mode'].values[0]
    
    # Compute the architecture. 
    arch = mut.thermo.SimpleRepression(R=R, ep_r=constants[g[1]], effector_conc=d['IPTGuM'],
                                       ka=Ka, ki=Ki, ep_ai=constants['ep_AI'], n_ns=constants['Nns'],
                                      n_sites=constants['n_sites']).bohr_parameter()
    if g[1] == 'O2':
        face = 'w'
    else:
        face = colors[g[0]]
    _ = collapse_ax.errorbar(arch, d['mean'], d['sem'], fmt='o', color=colors[g[0]], lw=1,
                            markerfacecolor=face, ms=4, linestyle='none')
    
# ###############################
# COLLAPSE CURVE
# ###############################
master_curve = (1 + np.exp(-bohr_range))**-1
collapse_ax.plot(bohr_range, master_curve, color='k', lw=1)

