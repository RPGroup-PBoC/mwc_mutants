# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import pandas as pd
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
colors = mut.viz.color_selector('pboc')
constants = mut.thermo.load_constants()
mut.viz.plotting_style()
    
# Load the sampling statistics
data = pd.read_csv('../../data/csv/epRA_sbc.csv')

# Compute the prior distribution functions
epRA_range = np.linspace(-40, 10, 500)
sigma_range = np.linspace(0, 0.5, 500)
epRA_pdf = scipy.stats.norm(-12, 6).pdf(epRA_range)
epRA_cdf = scipy.stats.norm(-12, 6).cdf(epRA_range)
sigma_pdf = scipy.stats.halfnorm(0, 0.1).pdf(sigma_range)
sigma_cdf = scipy.stats.halfnorm(0, 0.1).cdf(sigma_range)

# Define the bounds of the 95\% of the uniform distribution. 
n_bins = 20
n_sims = 10000
perc_low = scipy.stats.binom.ppf(0.005, 800, n_bins/ 400)
perc_high = scipy.stats.binom.ppf(0.095, 800, n_bins / 400)


# Instantiate the figure canvas
fig, ax = plt.subplots(2, 2, figsize=(7, 5))

# Apply special formatting. 
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=6.5)
    a.yaxis.set_tick_params(labelsize=6.5)

# Add labels. 
ax[0, 0].set_xlabel('DNA binding energy [$k_BT$]', fontsize=6.5)

ax[0, 1].set_xlabel('DNA binding energy [$k_BT$]', fontsize=6.5)
for i in range(2):
    ax[i, 0].set_ylabel(r'true $\Delta\varepsilon_{RA}$ [$k_BT$]', fontsize=6.5)
    ax[i, 0].set_xlabel(r'inferred $\Delta\varepsilon_{RA}$ [$k_BT$]', fontsize=6.5)
    ax[i, 1].set_ylabel('cumulative distribution', fontsize=6.5)

# # Plot the analytical prior distributions. 
ax[0, 0].plot(epRA_range, epRA_range, 'k-', lw=1.5, label=r'$\mathcal{N}(-12, 6)$')
ax[1, 0].plot(sigma_range, sigma_range, 'k-', lw=1.5, label=r'$\mathcal{N}(-12, 6)$')
# ax[0, 1].plot(epRA_range, epRA_cdf, 'k-', lw=1.5, label=r'$\mathcal{N}(-12, 6)$')
# ax[1, 0].plot(sigma_range, sigma_pdf, 'k-', lw=1.5, label=r'$\mathcal{HN}(0, 0.1)$')
# ax[1, 1].plot(sigma_range, sigma_cdf, 'k-', lw=1.5, label=r'$\mathcal{HN}(0, 0.1)$')

# Plot the simulated prior distributions. 
ax_ind = {'ep_RA': 0, 'sigma': 1}
n_bins = np.linspace(-40, 10, 80)
for g, d in data.groupby('param'):
    # Compute the histograms 
    ax[ax_ind[g], 0].plot(d['post_mean'], d['ground_truth'], '.', color=colors['red'], markersize=1)
    ax[ax_ind[g], 0].plot(d['post_mean'], d['ground_truth'], '.', color=colors['red'], markersize=1)
    # Compute the CDF
    x, y = np.sort(d['post_mean']), np.arange(0, len(d), 1) / len(d)
    _ = ax[ax_ind[g], 1].step(x, y, color=colors['blue'], lw=2)
    x, y = np.sort(d['ground_truth']), np.arange(0, len(d), 1) / len(d)
    _ = ax[ax_ind[g], 1].step(x, y, color=colors['red'], lw=2)    
    
plt.tight_layout()
plt.savefig('./test.pdf')
