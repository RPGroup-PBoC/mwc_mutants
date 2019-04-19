# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.viz
import mut.thermo
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load all of the data and samples
data = pd.read_csv('../../data/csv/pathological_F_data.csv')
samples = pd.read_csv('../../data/csv/pathological_F_samples.csv')
stats = pd.read_csv('../../data/csv/pathological_F_stats.csv')

# Define constants for plotting
bohr_range = np.linspace(-10, 10, 100)

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(2, 2)
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
ax[0, 0].set_xscale('log')
ax[1, 0].set_xscale('log')

# ##############################################################################
# TRUE DATA
# ##############################################################################
ax[0, 1].plot(bohr_range, (1 + np.exp(-bohr_range))**-1, 'k-', label='master curve')
ax[1, 1].plot(bohr_range, (1 + np.exp(-bohr_range))**-1, 'k-', 
             label='master curve')
ax[0, 0].plot(data['IPTGuM'].unique(), data['fc_mu'].unique(), 'k.',
             markerfacecolor='w', label='true fold-change')
ax[1, 0].plot(data['IPTGuM'].unique(), data['fc_mu'].unique(), 'k-', 
              label='true fold-change', zorder=1000)
ax[0, 1].plot(data['bohr'].unique(), data['fc_mu'].unique(), 'k.', 
              markerfacecolor='w', label='true free energy')

# ##############################################################################
# SIMULATED DATA 
# ##############################################################################
for g, d in data.groupby(['IPTGuM']):
    ax[1, 0].plot(d['IPTGuM'], d['fold_change'], '.', color=colors['red'], 
                alpha=0.5)

# ##############################################################################
# INFERRED DATA
# ##############################################################################
for g, d in stats.groupby('IPTGuM'):
    fc = d[d['parameter']=='fc_mu']
    bohr = d[d['parameter']=='empirical_bohr']
    ax[1, 1].plot(bohr['median'], fc['median'], '.', color=colors['blue'])
    ax[1, 0].plot(fc['IPTGuM'], fc['median'], '.', color=colors['blue'])

plt.tight_layout()    

# ##############################################################################
# FIGURE 2: DELTA F AND DELTA FC
# ##############################################################################
fig, ax = plt.subplots(2, 2, figsize=(6, 4))
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)

# ##############################################################################
# ZERO LINES
# ##############################################################################
ax[0, 0].plot(np.arange(0, 1, 0.1), np.arange(0, 1, 0.1), 'k-')

# ##############################################################################
# ERROR IN FC
# ##############################################################################
for g, d in stats.groupby('IPTGuM'):
    fc_mu = d[d['parameter']=='fc_mu']
    ax[0, 0].plot(fc_mu['true_mu'], fc_mu['median'], '.', color=colors['blue'])
    ax[0, 0].vlines(fc_mu['median'], fc_mu['hpd_min'], fc_mu['hpd_max'],  
                    color=colors['blue'])

ax[0, 0].set_xlim([0.9, 1.01])
ax[0,0].set_ylim([0.9, 1.01])

# ##############################################################################
# ERROR IN F
# ##############################################################################
for g, d in stats.groupby('IPTGuM'):
    dbohr = d[d['parameter']=='empirical_bohr']
    ax[0, 1].plot(dbohr['true_bohr'], dbohr['median'], '.', color=colors['blue'])

