# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.thermo
import mut.stats
import mut.viz
color = mut.viz.color_selector('mut')
pboc = mut.viz.color_selector('pboc')
constants = mut.thermo.load_constants()
mut.viz.plotting_style()

#  Load the data
data = pd.read_csv('../../data/csv/summarized_data.csv')
data = data[(data['class'] == 'DNA') & (data['operator']=='O2')]
samples = pd.read_csv('../../data/csv/Fig2_O2_DNA_binding_energy_samples.csv')
stats = pd.read_csv('../../data/csv/Fig2_O2_DNA_binding_energy_stats.csv')

# Determine the unique repressor copy numbers
reps = np.sort(data['repressors'].unique())
c_range = np.logspace(-2, 4, 200)
fig, ax = plt.subplots(len(reps), len(reps), figsize=(7,7), sharex=True, sharey=True) 

# Format the axes
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=8) 
    a.yaxis.set_tick_params(labelsize=8)
    a.set_xscale('log')
    a.set_xlim([1E-8, 1E-2])

# Add appropriate labels
for i in range(len(reps)):
    ax[i, 0].set_ylabel('fold-change', fontsize=8)
    ax[-1, i].set_xlabel('IPTG [M]', fontsize=8)
    ax[i, 0].text(-0.7, 0.5, '$R = $' + str(int(reps[i])), fontsize=8, backgroundcolor=pboc['pale_yellow'],
                 transform=ax[i,0].transAxes, rotation='vertical')
    ax[0, i].set_title('$R = $' + str(int(reps[i])), fontsize=8, backgroundcolor=pboc['pale_yellow'], y=1.08)
    ax[i, i].set_facecolor('#e4e7ec')
for i in range(4):
    ax[-1, i].set_xticks([1E-6, 1E-3])
    
# Add predictor titles
fig.text(-0.07, 0.55, 'predictor strain', fontsize=8, backgroundcolor='#E3DCD0', rotation='vertical')
fig.text(0.435, 0.95, 'comparison strain', fontsize=8, backgroundcolor='#E3DCD0', rotation='horizontal')
    
# Plot the data. 
for g, d in data.groupby(['mutant']):
    g = g.upper()
    for i, _ in enumerate(reps):
        for j, _ in enumerate(reps):
            _d = d[d['repressors'] == reps[j]]
            if i == j:
                face = 'w'
                edge = color[g]
            else:
                face = color[g]
                edge = color[g]
            _ = ax[i, j].errorbar(_d['IPTGuM'] / 1E6, _d['mean'], _d['sem'], markerfacecolor=face,
                                 markeredgecolor=edge, color=edge, lw=0.15, linestyle='none', fmt='o',
                                 ms=2.5, label=g)
           
            # Plot the best-fit lines. 
            for k, m in enumerate(data['mutant'].unique()):
                _d = data[(data['mutant']=='m') & (data['repressors'] == reps[j])]
                
                # Get the binding energies. 
                epRA = stats[stats['parameter']=='ep_RA.{}.{}'.format(m, int(reps[i]))].values[0][1:].astype(float)
                
                # Compute the fold-change
                epRA_mesh, c_mesh = np.meshgrid(epRA, c_range)
                fc = mut.thermo.SimpleRepression(R=reps[j], ep_r=epRA_mesh, ka=constants['Ka'], 
                                                 ki=constants['Ki'], ep_ai=constants['ep_AI'],
                                                effector_conc=c_mesh, n_sites=constants['n_sites'],
                                                n_ns=constants['Nns']).fold_change()
                
                # Plot the fit. 
                _ = ax[i, j].plot(c_range / 1E6, fc[:, 0], color=color[m], lw=0.75) 
                _ = ax[i, j].fill_between(c_range / 1E6, fc[:, 1], fc[:, 2], color=color[m], alpha=0.2) 

            
_  = ax[0, 3].legend(fontsize=8, bbox_to_anchor=(1.04, 0.95))
plt.subplots_adjust(wspace=0.05, hspace=0.05)
plt.savefig('FigS2_pairwise_predictions.pdf', bbox_inches='tight')