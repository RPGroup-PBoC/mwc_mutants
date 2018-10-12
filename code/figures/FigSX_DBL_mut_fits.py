# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.thermo
import mut.stats
import mut.viz
constants = mut.thermo.load_constants()
colors = mut.viz.color_selector('mut')
pboc = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load the summarized data and sampling statistics
data = pd.read_csv('../../data/csv/summarized_data.csv')
data = data[data['class']=='DBL']
samples = pd.read_csv('../../data/csv/Fig5_O2_DBL_samples.csv')
modes = samples.iloc[np.argmax(samples['logp'].values)]
c_range = np.logspace(-2, 4, 200)

# Set up the figure canvas
fig, ax = plt.subplots(3, 3, sharex=True, sharey=True, figsize=(5, 5))
for a in ax.ravel():
    a.set_xscale('log')
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    
rows = ['Y20I', 'Q21M', 'Q21A']
cols = ['Q294K', 'F164T', 'Q294V']

for i, dna in enumerate(rows):
    for j, ind in enumerate(cols):
        
        # Isolate the data and plot. 
        _data = data[data['mutant']==f'{dna}-{ind}']
        _ = ax[i, j].errorbar(_data['IPTGuM'] / 1E6, _data['mean'], _data['sem'], ms=3, lw=1, linestyle='none',
                         color=colors[f'{dna}-{ind}'], fmt='o')
        
        # Compute the line of best fit.
        m = f'{dna}-{ind}'
        arch = mut.thermo.SimpleRepression(R=data['repressors'].unique(), ep_r=modes[f'ep_RA.{m}'],
                                          ka=modes[f'Ka.{m}'], ki=modes[f'Ki.{m}'], ep_ai=modes[f'ep_AI.{m}'],
                                          n_ns=constants['Nns'], n_sites=constants['n_sites'],
                                          effector_conc=c_range).fold_change()
       
        # Plot the fits. 
        _ = ax[i, j].plot(c_range / 1E6, arch, lw=1, color=colors[m])
        
plt.subplots_adjust(hspace=0.1, wspace=0.1)