# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.thermo
import mut.viz
pboc = mut.viz.color_selector('pboc')
color = mut.viz.color_selector('mut')
mut.viz.plotting_style()

# Load the data. 
data = pd.read_csv('../../data/csv/summarized_data.csv')
data = data[data['class']=='DBL']

# Load the predictions
predictions = pd.read_csv('../../data/csv/Fig4_O2_double_samples.csv')

# Define the rows and columns
rows = {d:i for i, d in enumerate(predictions['DNA_mutant'].unique())}
cols = {d:i for i, d in enumerate(predictions['IND_mutant'].unique())}


# Set up the plot. 
fig, ax = plt.subplots(3, 3, figsize=(3.5, 3.5), sharex=True, sharey=True)
for a in ax.ravel():
    a.set_xscale('log')
    a.set_xlim([1E-8, 1E-2])
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_xticks([1E-7, 1E-5, 1E-3])
    
# Add the appropriate labels.
for i in range(3):
    ax[-1, i].set_xlabel('IPTG [M]', fontsize=8)
    ax[i, 0].set_ylabel('fold-change', fontsize=8) 
for d, i in rows.items():
    ax[i, 0].text(-0.75, 0.55, '{}'.format(d), fontsize=8, rotation='vertical', backgroundcolor=pboc['pale_yellow'],
                 transform=ax[i,0].transAxes)
for d, i in cols.items():
    ax[0, i].set_title('{}'.format(d), fontsize=8, y=1.08, backgroundcolor=pboc['pale_yellow'])
    
# Plot the predictions and credible regions. 
for d, i in rows.items():
    for ind, j in cols.items():
        # Plot the credible regions
        pred = predictions[(predictions['DNA_mutant'] == d) & (predictions['IND_mutant'] == ind)] 
        _ = ax[rows[d], cols[ind]].plot(pred['IPTGuM'] / 1E6, pred['mode'], lw=0.5, color=color['{}-{}'.format(d, ind)])
        _ = ax[rows[d], cols[ind]].fill_between(pred['IPTGuM'] / 1E6, pred['hpd_min'], pred['hpd_max'],
                                               color=color['{}-{}'.format(d, ind)], alpha=0.3)
        
        # Plot the data. 
        _data = data[data['mutant']=='{}-{}'.format(d, ind)]
        _ = ax[rows[d], cols[ind]].errorbar(_data['IPTGuM'] / 1E6, _data['mean'], _data['sem'],
                                           lw=1, color=color['{}-{}'.format(d, ind)], ms=3, fmt='o',
                                           linestyle='none')
 

# plt.tight_layout()
plt.subplots_adjust(hspace=0.05, wspace=0.05)
plt.savefig('Fig4_epAI_correction.svg', bbox_inches='tight')