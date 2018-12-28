# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolor
import matplotlib.gridspec as gridspec
import mut.viz
import mut.stats
mut.viz.plotting_style()
colors = mut.viz.color_selector('pboc')

# Load the data and DNA binding energy estimation samples
data = pd.read_csv('../../data/csv/compiled_data.csv')
samples = pd.read_csv('../../data/csv/KaKi_only_samples.csv')

# Isolate a particular mutant and repressor copy number for illustration.
data = data[data['operator']=='O2']
samples = samples[samples['operator']=='O2']

# Define The figure canvas. 
fig = plt.figure(figsize=(7, 6))
gs = gridspec.GridSpec(14, 13)

# Define the subplots for the corner plots
# Q294V
ax0 = fig.add_subplot(gs[0:2, 0:2])
ax1 = fig.add_subplot(gs[2:4, 0:2])
ax2 = fig.add_subplot(gs[2:4, 2:4])
ax3 = fig.add_subplot(gs[4:6,0:2])
ax4 = fig.add_subplot(gs[4:6,2:4])
ax5 = fig.add_subplot(gs[4:6,4:6])

# Q294K
ax7 = fig.add_subplot(gs[8:10, 0:2])
ax8 = fig.add_subplot(gs[10:12, 0:2])
ax9 = fig.add_subplot(gs[10:12, 2:4])
ax10 = fig.add_subplot(gs[12:15,0:2])
ax11 = fig.add_subplot(gs[12:15,2:4])
ax12 = fig.add_subplot(gs[12:15,4:6])

# Define the PPC axes
ax6 = fig.add_subplot(gs[:6, 8:])
ax13 = fig.add_subplot(gs[8:, 8:])

# Apply special formatting
axes = [ax0, ax1, ax2, ax3,
       ax4, ax5, ax6, ax7,
       ax8, ax9, ax10, ax11,
       ax12, ax13]
for a in axes:
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    
# Group corner plot axes.
axes = {'Q294V': [ax0, ax1, ax2, ax3, ax4, ax5, ax6],
       'Q294K': [ax7, ax8, ax9, ax10, ax11, ax12, ax13]}
for m, ax in axes.items():
    # Set the correctly displayed axes
    ax[0].set_xticklabels([]) 
    ax[0].set_yticks([])
    ax[1].set_xticklabels([])
    ax[2].set_xticklabels([])
    ax[2].set_yticks([])
    ax[4].set_yticklabels([])
    ax[5].set_yticklabels([])
    
    # Set the labeling. 
    ax[1].set_ylabel('$K_A$ [µM]', fontsize=8)
    ax[3].set_ylabel('$\sigma$', fontsize=8)
    ax[3].set_xlabel('$K_I$ [µM]', fontsize=8)
    ax[4].set_xlabel('$K_A$ [µM]', fontsize=8)
    ax[5].set_xlabel('$\sigma$', fontsize=8)

# # Add panel labels. 
ax0.text(-0.6, 1.1, '(A)', fontsize=8, transform=ax0.transAxes)
ax6.text(-0.4, 1.04, '(B)', fontsize=8, transform=ax6.transAxes)
ax7.text(-0.6, 1.1, '(C)', fontsize=8, transform=ax7.transAxes)
ax13.text(-0.4, 1.04, '(D)', fontsize=8, transform=ax13.transAxes)

# Add the glyphs
for m, ax in axes.items():
    # Isolate the samples and data
    _samples = samples[samples['mutant']==m]
    _samples.dropna(axis=1, inplace=True)
    _data = data[data['mutant']==m]
    reps = [key for key in _samples.keys() if 'y_rep' in key]   
    
    # Plot the joint distributions. 
    cmap = plt.cm.get_cmap('magma')
    normed = mplcolor.Normalize(vmin=np.min(_samples['lp__']), vmax=np.max(_samples['lp__']))
    _colors = [cmap(normed(v)) for v in _samples['lp__']]
    _ = ax[1].scatter(x=_samples['Ki'], y=_samples['Ka'], c=_colors, marker='.',
                     s=0.05, rasterized=True)
    _ = ax[3].scatter(_samples['Ki'], _samples['sigma'], marker='.', c=_colors, s=0.05, alpha=0.5, rasterized=True)
    _ = ax[4].scatter(_samples['Ka'], _samples['sigma'], marker='.', c=_colors, s=0.05, alpha=0.5, rasterized=True)
    
    # Plot the marginal histograms
    ka_hist, ka_bins = np.histogram(_samples['Ka'], bins=20)
    ki_hist, ki_bins = np.histogram(_samples['Ki'], bins=20)
    sigma_hist, sigma_bins = np.histogram(_samples['sigma'], bins=20)
    _ = ax[0].step(ki_bins[:-1], ki_hist, color=colors['red'])
    _ = ax[0].fill_between(ki_bins[:-1], ki_hist, color=colors['light_red'],
                          step='pre')
    _ = ax[2].step(ka_bins[:-1], ka_hist, color=colors['red'])
    _ = ax[2].fill_between(ka_bins[:-1], ka_hist, color=colors['light_red'],
                          step='pre')
    _ = ax[5].step(sigma_bins[:-1], sigma_hist, color=colors['red'])
    _ = ax[5].fill_between(sigma_bins[:-1], sigma_hist, 
                           color=colors['light_red'], step='pre')

    # Plot the ppc ecdfs.
    x, y = np.sort(_samples[reps]), np.arange(0, len(reps), 1) / len(reps)
    lower = np.percentile(x, 2.5, axis=0)
    upper = np.percentile(x, 97.5, axis=0)
    ax[-1].step(lower, y, '-', color=colors['red'], lw=1, label='__nolegend__')
    ax[-1].step(upper, y, 'k-', color=colors['red'], lw=1, label='__nolegend__')
    ax[-1].fill_betweenx(y, lower, upper, color=colors['light_red'])

    ax[-1].step([], [], '-', color=colors['red'], label='sample $\sim \mathcal{N}(\mu, \sigma)$') 
    
    # Plot the observed fold-change ecdf. 
    x, y = np.sort(_data['fold_change']), np.arange(0, len(_data), 1) / len(_data)
    ax[-1].step(x, y, 'k-', label='observed')
    ax[-1].legend(fontsize=8, handlelength=1)
plt.savefig('../../figures/FigSX_KaKi_ppc.svg', bbox_inches='tight')





