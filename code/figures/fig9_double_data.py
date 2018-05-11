# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import matplotlib.pyplot as plt
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
DNA_idx = pd.DataFrame(DNA_chains.iloc[np.argmax(DNA_chains['lnprobability'])]).reset_index()
DNA_idx.columns = ['parameter', 'mode']
IND_idx = pd.DataFrame(IND_chains.iloc[np.argmax(IND_chains['lnprobability'])]).reset_index()
IND_idx.columns = ['parameter', 'mode']
modes = pd.concat([DNA_idx, IND_idx])

# %%
c_range = np.logspace(-8, -2, 500)
# Set up the figure canvas.
fig, ax = plt.subplots(3, 3, figsize=(5, 5), sharex=True, sharey=True)
DNA_muts = ['Y20I', 'Q21A', 'Q21M']
IND_muts = ['F164T', 'Q294K', 'Q294V']

# Set up the axes locator.
axes = {}
for i, dna in enumerate(DNA_muts):
    for j, ind in enumerate(IND_muts):
        axes['{}-{}'.format(dna, ind)] = ax[i, j]

# Loop through the axes and plot the data theoretical predictions.
grouped = fc_df.groupby(['mutant'])
for g, d in grouped: 

    # Compute and plotthe theoretical prediction.
    if 'Q294K' not in g:
        # Get the stats.
        DNA_mut, IND_mut = g.split('-')
        ep_R = modes[modes['parameter']=='{}_eps_r'.format(DNA_mut)]['mode'].values[0]
        ka = modes[modes['parameter']=='{}_Ka'.format(DNA_mut)]['mode'].values[0]
        ki = modes[modes['parameter']=='{}_Ki'.format(DNA_mut)]['mode'].values[0]

        # Set up the architecture.
        arch = mut.thermo.SimpleRepression(R=260, ep_r=ep_R, ep_ai=4.5, ka=ka/1E6, ki=ki/1E6, effector_conc=c_range)
    axis = axes[g]
    _ = axis.errorbar(d['IPTGuM'] / 1E6, d['mean'], d['sem'], ms=5, fmt='.', lw=1, color=colors[g])

ax = ax.ravel()
for a in ax:
    a.set_xscale('log')
    a.tick_params(labelsize=8)
plt.tight_layout()



# %%
modes