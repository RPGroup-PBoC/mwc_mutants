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
fig, ax = plt.subplots(1, 2, figsize=(6, 2.5))
for a in ax:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)


# ##############################################################################
# SIMULATED FOLD-CHANGE
# ##############################################################################
for g, d in data.groupby(['draw']):
    ax[0].plot(d['bohr'], d['fold_change'], '.', color=colors['red'])

# ##############################################################################
# INFERRED FOLD-CHANGE
# ##############################################################################
for g, d in stats.groupby('true_bohr'):
    fc = d[d['parameter']=='fc_mu']
    ax[0].plot(fc['true_bohr'], fc['median'], '.', color=colors['blue'])

# ##############################################################################
# COLLAPSE CURVE
# ##############################################################################
ax[0].plot(bohr_range, (1 + np.exp(-bohr_range))**-1, 'k-')

# ##############################################################################
# DELTA BOHR
# ##############################################################################
for g, d in stats.groupby(['true_bohr']):
    delta_bohr = d[d['parameter']=='empirical_bohr']
    corr = d[d['parameter']=='delta_bohr']
    ax[1].plot(delta_bohr['true_bohr'], corr['median'], '.', color=colors['blue'])
    # ax[1].plot(corr['true_bohr'], corr['median'], '.', color=colors['green'])
    # ax[1].plot(delta_bohr['true_borh'], delta_bohr['median'], color=colors['blue'])


