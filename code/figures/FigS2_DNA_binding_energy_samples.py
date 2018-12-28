# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mut.viz
import mut.stats
mut.viz.plotting_style()
colors = mut.viz.color_selector('pboc')

# Load the data and DNA binding energy estimation samples
data = pd.read_csv('../../data/csv/compiled_data.csv')
samples = pd.read_csv('../../data/csv/DNA_binding_energy_samples.csv')

# Isolate a particular mutant and repressor copy number for illustration.
data = data[(data['mutant']=='Q21M') & (data['repressors']==260) & 
            (data['operator']=='O2')]

samples = samples[(samples['mutant']=='Q21M') & (samples['repressors']==260) &
                 (samples['operator']=='O2')]

# Drop unnecessary nans. 
samples.dropna(axis=1, inplace=True)

# Get the names of the draws. 
reps = [key for key in samples.keys() if 'y_rep' in key]


# Define The figure canvas. 
fig = plt.figure(figsize=(6, 3))
gs = gridspec.GridSpec(4, 9)
ax0 = fig.add_subplot(gs[0:2, 0:2])
ax1 = fig.add_subplot(gs[2:4, 0:2])
ax2 = fig.add_subplot(gs[2:4, 2:4])
ax3 = fig.add_subplot(gs[:,5:])

# Apply special formatting
axes = [ax0, ax1, ax2, ax3]
for a in axes:
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
ax0.set_yticks([])
ax2.set_yticks([])
ax0.set_xticklabels([])

# Add appropriate labes. 
ax1.set_xlabel(r'$\Delta\varepsilon_{RA}$ [$k_BT$]', fontsize=8)
ax1.set_ylabel(r'$\sigma$', fontsize=8)
ax2.set_xlabel(r'$\sigma$', fontsize=8)
ax3.set_xlabel('fold-change', fontsize=8)
ax3.set_ylabel('cumulative distribution', fontsize=8)

# Add panel labels. 
ax0.text(-0.4, 1.1, '(A)', fontsize=8, transform=ax0.transAxes)
ax3.text(-0.4, 1.04, '(B)', fontsize=8, transform=ax3.transAxes)
    
# Plot the joint distribution. 
_ = ax1.plot(samples['ep_RA'], samples['sigma'], ',', 
             color=colors['red'], rasterized=True, ms=0.3)

# Plot the marginal histograms. 
epRA_hist, epRA_bins = np.histogram(samples['ep_RA'], bins=20)
sigma_hist, sigma_bins = np.histogram(samples['sigma'], bins=20)
_ = ax0.step(epRA_bins[:-1], epRA_hist, color=colors['red'])
_ = ax0.fill_between(epRA_bins[:-1], epRA_hist, color=colors['light_red'],
                     step='pre')
_ = ax2.step(sigma_bins[:-1], sigma_hist, color=colors['red'])
_ = ax2.fill_between(sigma_bins[:-1], sigma_hist, color=colors['light_red'],
                     step='pre')

# Plot the ppc ecdfs.
x, y = np.sort(samples[reps]), [np.arange(0, len(reps), 1) / len(reps)] * len(samples)
ax3.step(x, y, '-', color=colors['red'], lw=0.5, alpha=0.8, label='__nolegend__')
ax3.step([], [], '-', color=colors['red'], label='sample $\sim \mathcal{N}(\mu, \sigma)$') 
    
# Plot the observed fold-change ecdf. 
x, y = np.sort(data['fold_change']), np.arange(0, len(data), 1) / len(data)
ax3.step(x, y, 'k-', lw=1, label='observed')
ax3.legend(fontsize=8, loc='upper left')
plt.savefig('../../figures/FigS2_epRA_ppc.pdf', bbox_inches='tight')
