#-*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
sys.path.insert(0, '../../')
import mut.viz
import mut.stats
import mut.thermo
constants = mut.thermo.load_constants()
colors = mut.viz.color_selector('pboc')
bright_colors = {c:i for c, i in colors.items() if '_' not in c}
mut.viz.plotting_style()

# Load the summarized data and samples
data = pd.read_csv('../../data/csv/summarized_data.csv')
data = data[(data['class']=='DNA') & (data['operator']=='O2')].copy()
global_stats = pd.read_csv('../../data/csv/FigS1_O2_DNA_binding_energy_global_stats.csv')
indiv_stats = pd.read_csv('../../data/csv/Fig2_O2_DNA_binding_energy_stats.csv')

# Isolate the most-likely parameter values
modes = global_samples.iloc[np.argmax(global_samples['logp'].values)]

# Define the concentration range for theory curves. 
c_range = np.logspace(-2, 4, 200)

# Set up the figure canvas
fig, ax = plt.subplots(2, 3, figsize=(5.75, 3.5))

# Define the mutant axes and repressor colors.
mut_ax = {'Y20I':ax[:, 0], 'Q21A': ax[:, 1], 'Q21M':ax[:, 2]}
rep_colors = {r:list(bright_colors.values())[i] for i, r in enumerate(data['repressors'].unique())}
rep_idx = {r:i for i, r in enumerate(data['repressors'].unique())}

# Format all axes labels
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
   
fig.text(0.01, 0.95, '(A)', fontsize=8)
fig.text(0.38, 0.95, '(B)', fontsize=8)
fig.text(0.63, 0.95, '(C)', fontsize=8)

for i in range(3):
    # Format the axes
    ax[1, i].set_xticks([0, 1, 2, 3])
    ax[1, i].set_xticklabels(np.sort(data['repressors'].unique().astype(int)), rotation='vertical')
    ax[0, i].set_xscale('log')
    ax[1, i].set_xlim([-0.5, 3.5])
    ax[0, i].set_ylim([-0.05, 1.2])
    ax[1, i].set_ylim([-9, -16.5])
    ax[0, i].set_xlim([1E-8, 1E-2])
    ax[0, i].set_xticks([1E-6, 1E-3])
    if i > 0:
        ax[0, i].set_yticklabels([])
        ax[1, i].set_yticklabels([])
    if i == 0: 
        # Set conditional labeling
        ax[1, i].set_ylabel('binding energy [$k_BT$]', fontsize=8)
        ax[0, i].set_ylabel('fold-change', fontsize=8) 
    
    # Add appropriate labels and titles
    ax[1, i].set_xlabel('rep. / cell', fontsize=8) 
    ax[0, i].set_xlabel('IPTG [M]', fontsize=8)
    ax[0, i].set_title(list(mut_ax.keys())[i], backgroundcolor=colors['pale_yellow'], fontsize=8,
                      y=1.02)
 
# Plot the data and global fits. 
for g, d in data.groupby(['mutant', 'repressors']):
    mut_ax[g[0]][0].errorbar(d['IPTGuM'] / 1E6, d['mean'], d['sem'], ms=2, linestyle='none',
                            lw=1, color=rep_colors[g[1]], label=int(g[1]), fmt='o')
    
    # Extract the stats
    epRA_stats = global_stats[global_stats['parameter']==f'ep_RA.{g[0]}'].values[0][1:].astype(float)
    epR_mesh, c_mesh = np.meshgrid(epRA_stats, c_range)
    
    # Compute the theoretical fold-change
    fc_theo = mut.thermo.SimpleRepression(R=g[1], ep_r=epR_mesh, effector_conc=c_mesh,
                                         ka=constants['Ka'], ki=constants['Ki'],
                                         n_sites=constants['n_sites'], n_ns=constants['Nns'],
                                         ep_ai=constants['ep_AI']).fold_change()
    
    # Plot the theory curves. 
    _ = mut_ax[g[0]][0].plot(c_range / 1E6, fc_theo[:, 0], lw=1, color=rep_colors[g[1]], label='__nolegend__')
    _ = mut_ax[g[0]][0].fill_between(c_range/1E6, fc_theo[:, 1], fc_theo[:, 2], color=rep_colors[g[1]], 
                                    alpha=0.4)
    
    # Find the individual stats and plot. 
    indiv = indiv_stats[indiv_stats['parameter']==f'ep_RA.{g[0]}.{int(g[1])}'].values[0][1:].astype(float)
    _ = mut_ax[g[0]][1].plot(rep_idx[g[1]], indiv[0], 'o', markerfacecolor='w', markeredgecolor=rep_colors[g[1]],
                            ms=3)
    _ = mut_ax[g[0]][1].vlines(rep_idx[g[1]], indiv[1], indiv[2], lw=1, color=rep_colors[g[1]], label='__nolegend__')

    
# Plot the global-fit value of the binding energy as a horizontal line
h_space = np.linspace(-0.5, 3.5, 200)
for m, a in mut_ax.items():
    stats = global_stats[global_stats['parameter']==f'ep_RA.{m}'].values[0][1:].astype(float)
    _ = a[1].hlines(stats[0], -0.5, 3.5, color=colors['dark_green'], lw=1, zorder=1)
    _ = a[1].fill_between(h_space, stats[1], stats[2], color=colors['dark_green'], alpha=0.4)

# Add the llegend for the third plot
leg = ax[0, 2].legend(loc='upper left', title='rep. / cell', fontsize=7,
                     handletextpad=0.14)
leg.get_title().set_fontsize(7)
plt.subplots_adjust(hspace=0.4)    

# Save!
plt.savefig('FigS3_global_dna.pdf', bbox_inches='tight')

