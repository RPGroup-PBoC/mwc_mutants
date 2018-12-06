# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import mut.stats
import mut.thermo
import mut.viz
import matplotlib.pyplot as plt
colors = mut.viz.color_selector('mut')
constants = mut.thermo.load_constants()
mut.viz.plotting_style()

# Load the data and samples
data = pd.read_csv('../../data/csv/compiled_data.csv')
samples = pd.read_csv('../../data/csv/DBL_Hill_samples.csv')



# Define the figure canvas. 
fig, ax = plt.subplots(3, 3, figsize=(4,4),sharex=True, sharey=True)

# Format the axes. 
DNA = ['Y20I', 'Q21M', 'Q21A']
IND = ['Q294K','F164T', 'Q294V']
for i, dna in enumerate(DNA):
    ax[0, i].set_title(dna, fontsize=8, backgroundcolor='#FFedc0')
    for j, ind in enumerate(IND):
#         ax[i, 0].text(-0.4, 0.5, ind, fontsize=8, rotation='vertical',
        ax[i, j].set_xscale('log')
        ax[i, j].set_xlim([1E-8, 1E-2])
        ax[i, j].set_ylim([-0.02, 1.2])
        ax[i, j].yaxis.set_tick_params(labelsize=8)
           
# ###############################
# INDUCTION DATA 
# ###############################
for i, dna in enumerate(DNA):
    for j, ind in enumerate(IND):
        # Restrict the data. 
        _data = data[data['mutant']==f'{dna}-{ind}'].groupby(['IPTGuM']).mean().reset_index()
        ax[i, j].plot(_data['IPTGuM'] / 1E6, _data['fold_change'],marker = 'o',
                         color=colors[f'{dna}-{ind}'], ms=3, lw=0)
       
# ##############################
# BEST FIT AND CREDIBLE REGIONS
# ##############################
c_range = np.logspace(-2, 4, 200)
for i, dna in enumerate(DNA):
    for j, ind in enumerate(IND):
        if f'{dna}-{ind}' != 'Y20I-Q294K':
            # Isolate the samples. 
            samps = samples[samples['mutant']==f'{dna}-{ind}']
            
            # Compute the credible regions
            cred_region = np.zeros((2, len(c_range)))
            for k, c in enumerate(c_range):
                arch = samps['a'] + samps['b'] * ((c / samps['K'])**samps['n'] / (1 + (c / samps['K'])**samps['n']))
                cred_region[:, k] = mut.stats.compute_hpd(arch, 0.95)
                
            # Plot the credible region. 
            ax[i, j].fill_between(c_range/1E6, cred_region[0, :], cred_region[1, :], color=colors[f'{dna}-{ind}'],
                                 alpha=0.4)
plt.subplots_adjust(hspace=0.05, wspace=0.05) 
        
        
        
        