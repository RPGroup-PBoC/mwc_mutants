# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sys.path.insert(0, '../../')
import mut.viz
import mut.stats
import mut.thermo
colors = mut.viz.pub_style()

# Load the WT fits.
wt_fit_stats = pd.read_csv('../../data/csv/WT_global_fit_parameters.csv')
modes = {}
for i, p in enumerate(wt_fit_stats['parameter'].unique()):
    modes[p] = wt_fit_stats[wt_fit_stats['parameter'] == p]['mode'].values[0]

# Define experimental constants.
ka_wt = modes['ka'] / 1E6
ki_wt = modes['ki'] / 1E6
n_ns = 4.6E6
ep_ai = 4.5
c_range = np.logspace(-8, -2, 500)

# Define the DNA mutants
DNA_MUTS = ['Q21A', 'Q21M', 'Y20I']

# Load the entire data set along with the sampler traces.
data = pd.read_csv('../../data/csv/compiled_data.csv')
stats = pd.read_csv('../../data/mcmc/DNA_O2_epR_fit_statistics.csv')
global_stats = pd.read_csv('../../data/mcmc/DNA_O2_global_fit_statistics.csv')
DNA = data[(data['class'] == 'DNA') | (data['class'] == 'WT')]
DNA = DNA[DNA['fold_change'] >= 0]

# Group the dataframe by mutant, IPTG, and repressor and compute statistics.
grouped = DNA.groupby(['mutant', 'repressors', 'IPTGuM']
                      ).apply(mut.stats.compute_mean_sem)
mean_sem_df = pd.DataFrame(grouped).reset_index()

# Regroup now only on mutant and repressors.
grouped = mean_sem_df.groupby(['mutant', 'repressors'])

# %%
# Set up the figure canvas.
fig, ax = plt.subplots(2, 3, figsize=(8, 6))
ax[0, 0].axis('off')
for i in range(3):
    ax[1, i].set_xscale('log')
    ax[1, i].set_ylim([-0.05, 1.2])
    ax[1, i].set_xlim([1E-8, 1E-2])

# Plot the leakiness fits.
ax[0, 1].set_xscale('log')
ax[0, 1].set_yscale('log')
ax[0, 1].set_xlabel('repressors per cell')
ax[0, 1].set_ylabel('fold-change')
ax[0, 1].set_xlim([1, 1E4])
# Define the axes for the mutants.
axes = {'Y20I': ax[1, 0], 'Q21A': ax[1, 1], 'Q21M': ax[1, 2]}

# Define the colors for the repressors.
reps = DNA['repressors'].unique()
muts = ['Q21M', 'Q21A', 'Y20I']
color_choices = ['red', 'green', 'blue', 'purple']
rep_colors = {i: colors[j] for i, j in zip(reps, color_choices)}
mut_colors = {i: j for i, j in zip(
    muts, sns.color_palette('viridis', n_colors=4))}

# Plot the leakiness fits.
rep_range = np.logspace(0, 4, 500)
for i, m in enumerate(muts):
    # Extract the parameters.
    epR_vals = stats[stats['parameter']
                     == 'ep_R_{}'.format(m)][['mode', 'hpd_min', 'hpd_max'
                                              ]].values
    # mesh togehter the values
    R, EpR = np.meshgrid(rep_range, epR_vals)

    # Instantiate the architecture and compute the fold-change.
    arch = mut.thermo.SimpleRepression(R, EpR, effector_conc=0, ep_ai=ep_ai,
                                       ka=ka_wt, ki=ki_wt)
    fc_theo = arch.fold_change()

    _ = ax[0, 1].plot(rep_range, fc_theo[0, :], color=mut_colors[m], label=m)
    _ = ax[0, 1].fill_between(rep_range, fc_theo[1, :],
                              fc_theo[2, :], color=mut_colors[m], alpha=0.3)

_ = ax[0, 1].legend()

# Plot the theoretical curves.
for i, m in enumerate(DNA_MUTS):
    # Extract the parameters.
    epR_vals = stats[stats['parameter']
                     == 'ep_R_{}'.format(m)][['mode', 'hpd_min', 'hpd_max'
                                              ]].values

    # mesh togehter the values
    R, C, EpR = np.meshgrid(reps, c_range, epR_vals)

    # Instantiate the architecture and compute the fold-change.
    arch = mut.thermo.SimpleRepression(R, EpR, effector_conc=C, ep_ai=ep_ai,
                                       ka=ka_wt, ki=ki_wt)
    fc_theo = arch.fold_change()

    # Iterate through and plot.
    for j, r in enumerate(reps):
        axis = axes[m]
        _ = axis.plot(c_range, fc_theo[:, j, 0],
                      lw=1, color=rep_colors[r], label='__nolegend__')
        _ = axis.fill_between(
            c_range, fc_theo[:, j, 1], fc_theo[:, j, 2], color=rep_colors[r], alpha=0.5)

# Iterate through each mutant and plot the data.
for g, d in grouped:
    if g[0] == 'wt':
        pass
    else:
        # Figure out the correct axis.
        axis = axes[g[0]]
        # Plot the data with the correct colors.
        axis.errorbar(d['IPTGuM'] / 1E6, d['mean'],
                      d['sem'], ms=2, color=rep_colors[g[1]], linestyle='none',
                      fmt='o', linewidth=1, label=int(g[1]))
        axis.set_title(g[0] + ' $\Delta\varepsilon_{RA} = %1f\, k_BT$' % stats[stats['parameter'] == 'ep_R_{}'.format(
            g[0])]['mode'].values[0], backgroundcolor=colors['pale_yellow'], y=1.08, fontsize=10)
        axis.set_xlabel('IPTG [M]', fontsize=10)
        axis.set_ylabel('fold-change', fontsize=10)
        if g[0] == 'Q21M':
            _ = axis.legend(loc='upper left')

        # Plot the leakiness.
        leakiness = d[d['IPTGuM'] == 0]
        ax[0, 1].errorbar(leakiness['repressors'], leakiness['mean'],
                          leakiness['sem'], ms=2, color=mut_colors[g[0]], label=g[0],
                          fmt='o', linewidth=1)
plt.tight_layout()
plt.savefig('../../figures/DNA_fits.pdf')
