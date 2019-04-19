# -*- coding: utf-8 -*-cd code/figures
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import mut.viz
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load the sbc data
sbc_data = pd.read_csv('../../data/csv/IND_sbc.csv')

# ##############################################################################
# FIGURE INSTANTIATION 
# ##############################################################################
fig, ax = plt.subplots(2, 2, figsize=(6, 4))
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)

# Assign axes
axes = {'KaKi_only': 0, 'KaKi_epAI':1}

# Add labels
for i in range(2):
    ax[0, i].set_ylim([-5, 5])
    ax[0, i].set_xlim([-0.1, 1.1])
    ax[0, i].set_xlabel('shrinkage', fontsize=8)
    ax[0, i].set_ylabel('z-score', fontsize=8)
    ax[1, i].set_xlabel('rank statistic', fontsize=8)
    ax[1, i].set_ylabel('cumulative distribution', fontsize=8)
    ax[1, i].set_xlim([0, 600])
    ax[i, 0].set_title('$K_A$ and $K_I$ modified', fontsize=8, 
                    backgroundcolor=colors['pale_yellow'], y=1.08)
    ax[i, 1].set_title(r'$K_A$, $K_I$, and $\Delta\varepsilon_{AI}$ modified',
                     fontsize=8, backgroundcolor=colors['pale_yellow'], y=1.08)

# Add panel labels
fig.text(0, 0.95, '(A)', fontsize=8)
fig.text(0.5, 0.95, '(B)', fontsize=8)

# ##############################################################################
# SHRINKAGE AND Z-SCORE
# ##############################################################################
legend = {'Ka': '$K_A$', 'Ki':'$K_I$', 'ep_AI':r'$\Delta\varepsilon_{AI}$'}
for g, d in sbc_data.groupby(['model']):
    _ax = ax[0, axes[g]]
    ka = d[d['param']=='Ka']
    ki = d[d['param']=='Ki']
    print(g, np.sum(ka['shrinkage'] > 0.9) / len(ka))
    if g != 'KaKi_only':
        ep = d[d['param']=='ep_AI']
        _ax.plot(ep['shrinkage'], ep['z_score'], '.', color='k', 
                label=legend['ep_AI'], ms=1)
    
    _ax.plot(ka['shrinkage'], ka['z_score'], '.', color=colors['red'], 
            label=legend['Ka'], ms=1)
    _ax.plot(ki['shrinkage'], ki['z_score'], '.', color=colors['blue'], 
            label=legend['Ki'], ms=1)

# ##############################################################################
# TRUE UNIFORM DISTRIBUTION
# ##############################################################################
n_sim = sbc_data.sim_idx.max()
L = np.arange(0, n_sim, 1)
R = sbc_data.rank_ndraws.unique()

# Envelope of cdf 99%
y = scipy.stats.randint.cdf(L, 0, R)
std = np.sqrt(y * (1 - y) / n_sim)
low_perc = np.concatenate((scipy.stats.norm.ppf(0.005, y[:-1], std[:-1]), (1.0, )))
high_perc = np.concatenate((scipy.stats.norm.ppf(0.995, y[:-1], std[:-1]), (1.0, )))
for i in range(2):
    ax[1, i].fill_between(L, low_perc, high_perc, color='slategray', alpha=0.4,
                        label='__nolegend__')

# ##############################################################################
# RANK DISTRIBUTION
# ##############################################################################
for g, d in sbc_data.groupby(['model']):
    _ax = ax[1, axes[g]]
    ka = d[d['param']=='Ka']
    ki = d[d['param']=='Ki']
    ka_x = np.sort(ka['rank'])
    ki_x = np.sort(ki['rank'])
    y = np.arange(0, len(ka), 1) / len(ka)
    if g != 'KaKi_only':
        ep = d[d['param']=='ep_AI']
        ep_x = np.sort(ep['rank'])
        _ax.step(ep_x, y, 'k', label=r'$\Delta\varepsilon_{AI}$')
    _ax.step(ka_x, y, color=colors['red'], label='$K_A$')
    _ax.step(ki_x, y, color=colors['blue'], label='$K_I$')

# ##############################################################################
# LEGEND DISTRIBUTION
# ##############################################################################
for a in ax.ravel():
    leg = a.legend(title='parameter', fontsize=6, markerscale=3, 
                   labelspacing=0.01, loc='upper left')
    leg.get_title().set_fontsize(6)
plt.tight_layout()
plt.savefig('../../figures/FigSX_IND_sbc.pdf')