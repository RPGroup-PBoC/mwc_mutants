# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.viz
import mut.thermo
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load all of the data and samples
data = pd.read_csv('../../data/csv/pathological_F_data.csv')
samples = pd.read_csv('../../data/csv/pathological_F_samples.csv')
stats = pd.read_csv('../../data/csv/pathological_F_stats.csv')

# Define constants for plotting
bohr_range = np.linspace(-10, 10, 100)

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(1, 2, figsize=(6, 3.5))
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)

ax[0].set_ylim([-0.2, 1.3])
ax[0].set_xlim([-8, 8])
ax[1].set_xlim([-8, 8])
ax[1].set_ylim([-8, 8])

# Add labels
ax[0].set_xlabel('free energy [$k_BT$]', fontsize=8)
ax[0].set_ylabel('fold-change', fontsize=8)
ax[1].set_xlabel('true free energy [$k_BT$]', fontsize=8)
ax[1].set_ylabel('inferred free energy [$k_BT$]', fontsize=8)

# Add panel labels
fig.text(0, 0.95, '(A)', fontsize=8)
fig.text(0.5, 0.95, '(B)', fontsize=8)
# ##############################################################################
# MASTER CURVE AND PERFECT AGREEMENT
# ##############################################################################
ax[0].plot(bohr_range, (1 + np.exp(-bohr_range))**-1, 'k-', lw=1, 
            label='master curve', zorder=1000)
ax[1].plot(bohr_range, bohr_range, 'k-', lw=1, label='perfect agreement',
            zorder=1000)

# ##############################################################################
#  SIMULATED DATA
# ##############################################################################
ax[0].plot(data['bohr'], data['fold_change'], '.', color=colors['red'], ms=1,
           alpha=0.5, label='simulated fold-change')

# ##############################################################################
# INFERRED DATA
# ##############################################################################
used_labels = []
_sig = []
for g, d in stats.groupby('true_bohr'):
    mu = d[d['parameter']=='fc_mu']
    sig = d[d['parameter']=='fc_sigma']
    _sig.append(sig['median'].values[0])
    bohr = d[d['parameter']=='empirical_bohr']

    if (mu['median'].values[0] < sig['median'].values[0]):
        color='tomato'
        label = 'µ < σ'
    elif ((1 - mu['median'].values[0]) < sig['median'].values[0]):
        color = 'rebeccapurple'
        label = '1 - µ < σ' 
    else:
        label = 'µ > σ ; 1 - µ > σ'
        color = colors['blue']  
    if label in used_labels:
        label = '__nolegend__'
    else:
        used_labels.append(label)
    if g != -8: 
        _label = '__nolegend__'
    else:
        _label = 'inferred µ'

    # Plot the inferred fold-change as a function of the bohr parameter
    ax[0].plot(g, mu['median'], 'o', color=colors['blue'], ms=2, 
                label=_label)
    ax[0].vlines(g, mu['hpd_min'], mu['hpd_max'], lw=1, color=colors['blue'],
                label='__nolegend__')
    
    # Plot the inferred bohr parameter as a function of the true value
    ax[1].plot(g, bohr['median'], '.', ms=3, color=color, 
                alpha=0.5, label=label)
    ax[1].vlines(g, bohr['hpd_min'], bohr['hpd_max'], lw=0.5, 
                color=color, label='__nolegend__')

# ##############################################################################
# PLOT THE SIGMA THRESHOLDS
# ##############################################################################
thresh = -np.log(np.mean(_sig)**-1 - 1)
ax[1].vlines(thresh, -8, 8, linestyle='--', color=colors['green'])
ax[1].vlines(-thresh, -8, 8, linestyle='--', color=colors['green'])

# ##############################################################################
#  ADD F SIGMA LABELS
# ##############################################################################
ax[1].text(3, -5.5, r'$-\log\left(\frac{1}{1 - \sigma} - 1\right)$', fontsize=8, 
            color=colors['green'])
ax[1].text(-6.5, 3, r'$-\log\left(\frac{1}{\sigma} - 1\right)$', fontsize=8, 
            color=colors['green'])

# ##############################################################################
# LEGENDS 
# ##############################################################################
ax[0].legend(loc='upper left', fontsize=6)
leg = ax[1].legend(title='inferred free energy', fontsize=6, handlelength=1)
leg.get_title().set_fontsize(6)
plt.tight_layout()
plt.savefig('../../figures/FigSX_delta_F_pathology.pdf', bbox_inches='tight')

