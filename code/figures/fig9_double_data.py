# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
sys.path.insert(0, '../../')
import mut.viz
import mut.stats
import mut.thermo
mut.viz.plotting_style()
colors = mut.viz.color_selector('mut')
pboc_colors = mut.viz.color_selector('pboc')

# Load the data
data = pd.read_csv('../../data/csv/compiled_data.csv')
single_muts = data[(data['class'] != 'WT') & (data['class'] != 'DBL')]

# Restrict only to double mutants.
dbl = data[data['class'] == 'DBL'].copy()

# Generate a mean fold-change dataframe for easy plotting.
grouped = dbl.groupby(['mutant', 'IPTGuM']).apply(mut.stats.compute_mean_sem)
fc_df = pd.DataFrame(grouped).reset_index()

# Load the chains from the single mutant fits.
DNA_chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_DNA_strict.csv')
IND_chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_IND_strict.csv')

# Compute the stats for each.
DNA_idx = pd.DataFrame(DNA_chains.iloc[np.argmax(
    DNA_chains['lnprobability'])]).reset_index()
DNA_idx.columns = ['parameter', 'mode']
IND_idx = pd.DataFrame(IND_chains.iloc[np.argmax(
    IND_chains['lnprobability'])]).reset_index()
IND_idx.columns = ['parameter', 'mode']
modes = pd.concat([DNA_idx, IND_idx])

# Set some constants.
c_range = np.logspace(-8, -2, 500)
ka_wt = 139E-6
ki_wt = 0.53E-6
ep_R_wt = -13.9

# %% Set up the somewhat complicated figure.
fig = plt.figure(figsize=(6.5, 3))

# Set up the axes.
gs = gridspec.GridSpec(3, 8)
axes = [[], [], []]
for i in range(3):
    for j in range(3):
        _ax = fig.add_subplot(gs[i, j])
        _ax.tick_params(labelsize=6)
        _ax.set_xscale('log')
        _ax.set_xlim([1E-8, 1E-2])
        _ax.set_ylim([-0.1, 1.1])
        _ax.xaxis.set_ticks([1E-6, 1E-3])
        _ax.xaxis.set_tick_params()
        _ax.yaxis.set_ticks([0, 0.5, 1])
        axes[i].append(_ax)
ax_bohr = fig.add_subplot(gs[:, 4:])

# Properly label and format the axes.
dna_muts = ['Y20I', 'Q21A', 'Q21M']
for i in range(3):
    axes[0][i].set_title(
        dna_muts[i], backgroundcolor=pboc_colors['pale_yellow'], fontsize=8, y=1.04)
    axes[0][i].set_xticklabels([])
    axes[1][i].set_xticklabels([])
    if i > 0:
        axes[0][i].set_yticklabels([])
        axes[1][i].set_xticklabels([])
        axes[1][i].set_yticklabels([])
        axes[2][i].set_yticklabels([])

ind_muts = ['F164T', 'Q294K', 'Q294V']
for i in range(3):
    axes[i][0].text(-1.1, 0.65, ind_muts[i], fontsize=8, backgroundcolor=pboc_colors['pale_yellow'],
                    transform=axes[i][0].transAxes, rotation='vertical')
    axes[i][0].set_ylabel('fold-change', fontsize=8)
    axes[-1][i].set_xlabel('IPTG [M]', fontsize=8)
ax_bohr.tick_params(labelsize=8)
ax_bohr.set_xlabel('Bohr parameter [$k_BT$]', fontsize=8)
ax_bohr.set_ylabel('fold-change', fontsize=8)
ax_bohr.set_xlim([-10, 10])

# Add panel labels.
fig.text(0, 1, '(A)', fontsize=8)
fig.text(0.45, 1, '(B)', fontsize=8)

# plot the double mutant data and predictions.
# Set up the axes locator.
mut_axes = {}
for i, dna in enumerate(dna_muts):
    for j, ind in enumerate(ind_muts):
        mut_axes['{}-{}'.format(dna, ind)] = axes[j][i]

# Loop through the axes and plot the data theoretical predictions.
grouped = fc_df.groupby(['mutant'])
for g, d in grouped:
    axis = mut_axes[g]

    # oompute and plotthe theoretical prediction.
    if 'Q294K' not in g:
        # Get the stats.
        DNA_mut, IND_mut = g.split('-')
        ep_R = modes[modes['parameter'] ==
                     '{}_eps_r'.format(DNA_mut)]['mode'].values[0]
        ka = modes[modes['parameter'] ==
                   '{}_Ka'.format(IND_mut)]['mode'].values[0]
        ki = modes[modes['parameter'] ==
                   '{}_Ki'.format(IND_mut)]['mode'].values[0]

        # Set up the architecture.
        arch = mut.thermo.SimpleRepression(
            R=260, ep_r=ep_R, ep_ai=4.5, ka=ka/1E6, ki=ki/1E6, effector_conc=c_range)
        fc_theo = arch.fold_change()
        _ = axis.plot(c_range, fc_theo, lw=0.75, color=colors[g])

    _ = axis.errorbar(d['IPTGuM'] / 1E6, d['mean'], d['sem'],
                      ms=2, fmt='.', lw=1, color=colors[g])

# Plot the collapse curve.
# All data.
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[(data['mutant'] != 'Q294K') & (data['mutant'] != 'Q294R')].copy()
grouped = data.groupby(['mutant', 'repressors', 'IPTGuM'])
bohr = np.linspace(-10, 10, 500)
theo = (1 + np.exp(-bohr))**-1
ax_bohr.plot(bohr, theo, 'k-', label='prediction', zorder=1000)

labeled = []
for g, d in grouped:
    split = g[0].split('-')
    if 'Q294K' not in split:
        if len(split) == 2:
            ep_R = modes[modes['parameter'] ==
                         '{}_eps_r'.format(split[0])]['mode'].values[0]
            ka = modes[modes['parameter'] == '{}_Ka'.format(
                split[1])]['mode'].values[0] * 1E-6
            ki = modes[modes['parameter'] == '{}_Ki'.format(
                split[1])]['mode'].values[0] * 1E-6
        else:
            if split[0] in dna_muts:
                ep_R = modes[modes['parameter'] ==
                             '{}_eps_r'.format(split[0])]['mode'].values[0]
                ka = ka_wt
                ki = ki_wt
            else:
                ep_R = ep_R_wt
                ka = modes[modes['parameter'] == '{}_Ka'.format(
                    split[0])]['mode'].values[0] * 1E-6
                ki = modes[modes['parameter'] == '{}_Ki'.format(
                    split[0])]['mode'].values[0] * 1E-6

        # Set up the architecture.
        arch = mut.thermo.SimpleRepression(R=g[1], ep_r=ep_R, ka=np.abs(
            ka), ki=np.abs(ki), ep_ai=4.5, effector_conc=g[2]/1E6)
        bohr_param = arch.bohr_parameter()

        # Compute the mean and SEM of the fold_change.
        mean_fc = d['fold_change'].mean()
        sem_fc = d['fold_change'].std() / np.sqrt(len(d))
        if split[0] == 'wt':
            face = 'w'
            edge = 'k'
        else:
            face = colors[g[0]]
            edge = colors[g[0]]
        if g[0] not in labeled:
            label = g[0].upper()
            labeled.append(g[0])
        else:
            label = '__nolegend__'
        _ = ax_bohr.errorbar(bohr_param, mean_fc, sem_fc, fmt='.', color=colors[g[0].upper(
        )], lw=1, linestyle='none', alpha=0.8, markerfacecolor=face, markeredgecolor=edge, markeredgewidth=1, markersize=7, label=label)

ax_bohr.legend(bbox_to_anchor=(1.0, 1), fontsize=8)
plt.savefig('../../figures/doubles_data_collapse.pdf', bbox_inches='tight')
