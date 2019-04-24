# -*- coding; utf-8 -*- 
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
import seaborn as sns
constants = mut.thermo.load_constants()
colors = mut.viz.color_selector('pboc')
_colors = sns.color_palette('deep', n_colors=3)
mut.viz.plotting_style()

# Load and restrict the various data sets
data = pd.read_csv('../../data/csv/summarized_data.csv')
data = data[data['class']=='IND'].copy()
stats = pd.read_csv('../../data/csv/epAI_only_summary.csv')
stats = stats[stats['operator']=='O2'].copy()
bohr = pd.read_csv('../../data/csv/empirical_F_statistics.csv')
bohr = bohr[bohr['class']=='IND'].copy()

# Define constants for plotting 
c_range = np.logspace(-3, 4, 200)
c_range[0] = 0
bohr_range = np.linspace(-8, 8, 200)
F = (1 + np.exp(-bohr_range))**-1

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(3, 4, figsize=(6, 4))
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)

# Add appropriate scaling
for i in range(4):
    ax[0, i].set_xscale('symlog', linthreshx=0.006)
    ax[-1, i].set_xscale('symlog', linthreshx=0.006)
    ax[-1, i].set_ylim([-8, 8])
    ax[0, i].set_ylim([-0.2, 1.2])
    ax[1, i].set_ylim([-0.2, 1.2])
    ax[0, i].set_xticks([0, 1E-2, 1E0, 1E2, 1E4]) 
    ax[0, i].set_yticks([0, 0.5, 1])
    ax[1, i].set_yticks([0, 0.5, 1])
    ax[-1, i].set_xticks([0, 1E-2, 1E0, 1E2, 1E4]) 
    ax[1, i].set_xlim([-8, 8])

# Define the axes
axes = {'F164T':0, 'Q294V':1, 'Q294K':2, 'Q294R':3}
op_colors = {'O1':_colors[0], 'O2':_colors[1], 'O3':_colors[2]}

# ##############################################################################
#  COLLAPSE CURVE
# ##############################################################################
for i in range(4):
    ax[1, i].plot(bohr_range, F, 'k-', lw=1)

# ##############################################################################
# FOLD CHANGE CURVES 
# ##############################################################################
for i, o in enumerate(('O1', 'O2', 'O3')):
    ep_r = constants[o]
    for m, a in axes.items():
        _stats = stats[(stats['mutant']==m) & (stats['parameter']=='ep_AI')][
                    ['hpd_min', 'hpd_max']].values
        _c, _ep = np.meshgrid(c_range,_stats) 
        arch = mut.thermo.SimpleRepression(R=260, ep_r=ep_r, ka=constants['Ka'],
                                            ki=constants['Ki'], ep_ai=_ep,
                                            effector_conc=_c).fold_change()
        ax[0, a].fill_between(c_range, arch[0, :], arch[1, :], color=op_colors[o],
        alpha=0.4)

# ##############################################################################
# FOLD-CHANGE DATA 
# ##############################################################################
for g, d in data.groupby(['mutant', 'operator']):
    _ax = ax[0, axes[g[0]]]
    if g[1] == 'O2':
        face = 'w'
    else:
        face = op_colors[g[1]]
    _ax.errorbar(d['IPTGuM'], d['mean'], d['sem'], fmt='.', markersize=5, 
        markerfacecolor=face, linestyle='none', color=op_colors[g[1]], capsize=1)

# ##############################################################################
# DELTA F DATA
# ##############################################################################
for g, d in bohr.groupby(['mutant', 'operator', 'IPTGuM']):
    _ax = ax[2, axes[g[0]]]
    _param = d[d['parameter']=='delta_bohr_corrected2']
    mu = d[d['parameter']=='fc_mu']['median'].values[0]
    sig = d[d['parameter']=='fc_sigma']['median'].values[0]
    if (mu < sig) | (1 - mu < sig):
        color = 'slategray'
        alpha = 0.25
    else:
        color = op_colors[g[1]]
        alpha = 1 
    if g[1] == 'O2':
        face = 'w'
    else:
        face = color

    _ax.plot(_param['IPTGuM'], _param['median'], '.', color=color, 
            markerfacecolor=face , alpha=alpha)
plt.tight_layout()

