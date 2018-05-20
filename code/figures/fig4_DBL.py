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
mut.viz.plotting_style()
colors = mut.viz.color_selector('mut')
pboc_colors = mut.viz.color_selector('pboc')

# Load the data.
data = pd.read_csv('../../data/csv/compiled_data.csv')
DNA_muts = ['Y20I', 'Q21A', 'Q21M']
IND_muts = ['F164T', 'Q294K', 'Q294V']

DNA_chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_DNA_strict.csv')
IND_chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_IND_strict.csv')

c_range = np.logspace(-8, -2, 500)
# %%
# Set up the figure canvas
fig, ax = plt.subplots(3, 3, figsize=(3.5, 3.5), sharex=True, sharey=True)
axes = ax.ravel()
for a in axes:
    a.set_xscale('log')
    a.xaxis.set_ticks([1E-7,  1E-5,  1E-3])
    a.set_xlim([1E-8, 1E-2])
    a.set_ylim([-0.15, 1.15])
    a.tick_params(labelsize=8)


for i in range(3):
    axes[i].set_title(DNA_muts[i], fontsize=8, backgroundcolor=pboc_colors['pale_yellow'], y=1.05)
    ax[i, 0].text(-0.7, 0.65, IND_muts[i], fontsize=8, backgroundcolor=pboc_colors['pale_yellow'], rotation='vertical',
    transform=ax[i, 0].transAxes)

ax[1, 0].set_ylabel('fold-change', fontsize=8)
ax[-1, 1].set_xlabel('IPTG (M)', fontsize=8)

for i, dna in enumerate(DNA_muts):
    for j, ind in enumerate(IND_muts):
        # Isolate the chains and the data for the double mutant.
        dbl = data[data['mutant']=='{}-{}'.format(dna, ind)]
        if ind != 'Q294K':
            dna_chain = DNA_chains['{}_eps_r'.format(dna)]
            ka_chain = IND_chains['{}_Ka'.format(ind)]
            ki_chain = IND_chains['{}_Ki'.format(ind)]
            dna_idx = np.argmax(DNA_chains['lnprobability'])
            ind_idx = np.argmax(IND_chains['lnprobability'])
            ep_r = dna_chain.values[dna_idx]
            ka = ka_chain.values[ind_idx]
            ki = ki_chain.values[ind_idx]

            # Compute the architecture. 
            arch = mut.thermo.SimpleRepression(R=260, ep_r=ep_r, ep_ai=4.5, ka=ka/1E6,
            ki=ki/1E6, effector_conc=c_range).fold_change()
            _ = ax[j, i].plot(c_range, arch, color=colors['{}-{}'.format(dna, ind)], lw=1)
        
        # plot the data.
        grouped = dbl.groupby(['IPTGuM'])
        for g, d in grouped:
            _ = ax[j, i].errorbar(g/1E6, d['fold_change'].mean(), d['fold_change'].std() / np.sqrt(len(d)),
            color=colors['{}-{}'.format(dna, ind)], lw=1, fmt='.', markersize=4)

plt.subplots_adjust(hspace=0.05, wspace=0.05)
plt.savefig('../../figures/double_plots.svg', bbox_inches='tight')
#%%
data