# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.viz
import mut.thermo
constants = mut.thermo.load_constants()
mut_colors = mut.viz.color_selector('mut')
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load the statistics
stats = pd.read_csv('../../data/csv/DNA_binding_energy_summary.csv')
data = pd.read_csv('../../data/csv/summarized_data.csv')
data = data[data['class']=='DNA'].copy()

# Define some constants
c_range = np.logspace(-2, 4, 200)
c_range[0] = 0
repressors = [60, 124, 260, 1220]
muts = ['Y20I', 'Q21A', 'Q21M']

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(4, 4, figsize=(5, 5), sharex=True, sharey=True)
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_xscale('symlog')

for i, r in enumerate(repressors):
    ax[0, i].set_title(f'$R =$ {int(r)}', fontsize=8, 
            backgroundcolor=colors['pale_yellow'], y=1.08)
    ax[i, 0].text(-0.7, 0.65, f'R={int(r)}', rotation='vertical', fontsize=8,
                  backgroundcolor=colors['pale_yellow'], 
                  transform=ax[i,0].transAxes)
    ax[i, 0].set_yticks([0, 0.5, 1])
    ax[i,0].set_ylim([-0.1, 1.1])
    ax[i, 0].set_ylabel('fold-change', fontsize=8)
    ax[-1, i].set_xlabel('IPTG [µM]', fontsize=8)

# Change the background color of the diagonals
for i in range(4):
    ax[i, i].set_facecolor('#D1D3D4')

# ##############################################################################
#  FITS
# ##############################################################################
for i, m in enumerate(muts):
    for j, r1 in enumerate(repressors):
        fit_strain = stats[(stats['repressors']==r) & (stats['mutant']==m)\
                            & (stats['parameter']=='ep_RA')]
        for k, r2 in enumerate(repressors):
            # Compute the fit
            ep, c = np.meshgrid(fit_strain[['hpd_min', 'hpd_max']].values,
                                 c_range)
            arch = mut.thermo.SimpleRepression(R=r2, ep_r=ep, ka=constants['Ka'],
                                              ki=constants['Ki'], 
                                              n_sites=constants['n_sites'],
                                              ep_ai=constants['ep_AI'], 
                                              effector_conc=c).fold_change()
            ax[j, k].fill_between(c_range, arch[:, 0], arch[:, 1], 
                                  color=mut_colors[m], alpha=0.5)

# ##############################################################################
# EXPERIMENTAL DATA
# ##############################################################################
mut_colors = {'Y20I': colors['blue'], 'Q21A':colors['purple'], 'Q21M':colors['red']}
for i, r1 in enumerate(repressors):
    for j, r2 in enumerate(repressors):
        for k, m in enumerate(muts):
            if i == j:
                face = 'w'
            else:
                face = mut_colors[m]

            # Isolate the correct mutant and repressor copy number
            _data = data[(data['mutant']==m) & (data['repressors']==r2)]
            ax[i, j].errorbar(_data['IPTGuM'], _data['mean'], _data['sem'], 
                            linestyle='none', fmt='.', ms=4, capsize=1, 
                            linewidth=1, markerfacecolor=face, color=mut_colors[m],
                            label=m)

plt.tight_layout()

