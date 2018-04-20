# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.viz
import mut.stats
import mut.thermo
colors = mut.viz.pub_style()

# Define experimental constants.
ka_wt = 139E-6
ki_wt = 0.53E-6
n_ns = 4.6E6
ep_ai = 4.5
c_range = np.logspace(-8, -2, 500)

# Define the DNA mutants
DNA_MUTS = ['Q21A', 'Q21M', 'Y20I']

# Load the entire data set along with the sampler traces.
data = pd.read_csv('../../data/csv/compiled_data.csv')
stats = pd.read_csv('../../data/mcmc/DNA_O2_epR_fit_statistics.csv')
DNA = data[(data['class'] == 'DNA') | (data['class'] == 'WT')]


def compute_mean_sem(df):
    """
    Computes the mean and standard error of the fold-change given a
    grouped pandas Series.
    """
    # Compute the properties
    mean_fc = df['fold_change'].mean()
    sem_fc = df['fold_change'].std() / np.sqrt(len(df))

    # Assemble the new pandas series and return.
    samp_dict = {'mean': mean_fc, 'sem': sem_fc}
    return pd.Series(samp_dict)


# Group the dataframe by mutant, IPTG, and repressor and compute statistics.
grouped = DNA.groupby(['mutant', 'repressors', 'IPTGuM']
                      ).apply(compute_mean_sem)
mean_sem_df = pd.DataFrame(grouped).reset_index()

# Regroup now only on mutant and repressors.
grouped = mean_sem_df.groupby(['mutant', 'repressors'])

# %%
# Set up the figure canvas.
fig, ax = plt.subplots(2, 3, figsize=(6, 4))
for i in range(3):
    ax[0, i].set_xscale('log')
    ax[0, i].set_ylim([-0.05, 1.2])
    ax[0, i].set_xlim([1E-8, 1E-2])

# Define the axes for the mutants.
axes = {'Y20I': ax[0, 0], 'Q21A': ax[0, 1], 'Q21M': ax[0, 2]}

# Define the colors for the repressors.
reps = DNA.repressors.unique()
color_choices = ['red', 'green', 'blue', 'purple']
rep_colors = {i: colors[j] for i, j in zip(reps, color_choices)}

# Plot the theoretical curves.
for i, m in enumerate(DNA_MUTS):
    # Extract the parameters.
    epR_vals = stats[stats['parameter']
                     == 'ep_R_{}'.format(m)][['mode', 'hpd_min', 'hpd_max'
                     ]].values

    # mesh togehter the values
    R, C, EpR= np.meshgrid(reps, c_range, epR_vals)

    # Instantiate the architecture and compute the fold-change.
    arch= mut.thermo.SimpleRepression(R, EpR, effector_conc=C, ep_ai=ep_ai,
    ka=ka_wt, ki=ki_wt)
    fc_theo= arch.fold_change()

    # Iterate through and plot.
    for j, r in enumerate(reps):
        axis=axes[m]
        _ = axis.plot(c_range, fc_theo[:, j, 0], lw=1.5, color=rep_colors[r], label=int(r))
        _ = axis.fill_between(c_range, fc_theo[:, j, 1], fc_theo[:, j, 2], color=rep_colors[r], alpha=0.5)
# Iterate through each mutant and plot the data.
for g, d in grouped:
    if g[0] == 'wt':
        pass
    else:
        # Figure out the correct axis.
        axis= axes[g[0]]
        # Plot the data with the correct colors.
        axis.errorbar(d['IPTGuM'] / 1E6, d['mean'],
                      d['sem'], ms=3, color=rep_colors[g[1]], linestyle='none',
                      fmt='o')


plt.tight_layout()
