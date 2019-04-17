# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import seaborn as sns
import mut.viz
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load the posterior check data. 
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[(data['mutant']=='Q294K') & (data['operator']=='O2')]
kaki_data = pd.read_csv('../../data/csv/IND_posterior_predictive_samples.csv')
kaki_data = kaki_data[kaki_data['model']=='KaKi_only']
kaki_samples = pd.read_csv('../../data/csv/KaKi_only_samples.csv')
kaki_samples = kaki_samples[(kaki_samples['mutant']=='Q294K') &
                             (kaki_samples['operator']=='O2')]
cmap = sns.color_palette('viridis', n_colors=len(kaki_samples))
kaki_samples.sort_values('lp__', inplace=True)
kaki_samples['color'] = cmap

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig = plt.figure(figsize=(7, 3))
gs = gridspec.GridSpec(6, 11)
ax1 = fig.add_subplot(gs[0:2, 0:2])
ax2 = fig.add_subplot(gs[2:4, 0:2])
ax3 = fig.add_subplot(gs[2:4, 2:4])
ax4 = fig.add_subplot(gs[4:6, 0:2])
ax5 = fig.add_subplot(gs[4:6, 2:4])
ax6 = fig.add_subplot(gs[4:6, 4:6])
ax7 = fig.add_subplot(gs [:, 7:])
ax = [ax1, ax2, ax3, ax4, ax5, ax6, ax7]

# Apply special formatting where needed
for a in ax:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
for a in [ax1, ax2, ax3]:
    a.set_xticklabels([])
for a in [ax1, ax3, ax5, ax6]:
    a.set_yticks([])

# Adjust scaling
ax7.set_xscale('symlog', linthreshx=1E-2)

# Add labels
ax4.set_xlabel('$K_A$ [µM]', fontsize=8)
ax2.set_ylabel('$\sigma$ ', fontsize=8)
ax4.set_ylabel('$K_I$ [µM]', fontsize=8)
ax5.set_xlabel('$\sigma$', fontsize=8)
ax6.set_xlabel('$K_I$ [µM]', fontsize=8)
ax7.set_xlabel('IPTG [µM]', fontsize=8)
ax7.set_ylabel('fold-change', fontsize=8)

# Add panel labels
fig.text(0.02, 0.95, '(A)', fontsize=8)
fig.text(0.5, 0.95, '(B)', fontsize=8)
# ##############################################################################
# JOINT DISTRIBUTIONS
# ##############################################################################
ax2.scatter(kaki_samples['Ka'].values, kaki_samples['sigma'].values, marker=',',
            c=kaki_samples['color'].values,  s=0.1)
ax4.scatter(kaki_samples['Ka'].values, kaki_samples['Ki'].values, marker=',',
            c=kaki_samples['color'].values, s=0.1)
ax5.scatter(kaki_samples['sigma'].values, kaki_samples['Ki'].values, marker=',',
            c=kaki_samples['color'].values, s=0.1)
            
# ##############################################################################
# MARGINAL DISTRIBUTIONS
# ##############################################################################
ax1.hist(kaki_samples['Ka'], bins=30, histtype='stepfilled', color=cmap[0], 
         alpha=0.5)
ax3.hist(kaki_samples['sigma'], bins=30, histtype='stepfilled', color=cmap[0], 
         alpha=0.5)
ax6.hist(kaki_samples['Ki'], bins=30, histtype='stepfilled', color=cmap[0],
         alpha=0.5)

# ##############################################################################
# PERCENTILES
# ##############################################################################
percs = [99, 95, 80, 50, 20, 10, 5]
perc_cmap = sns.color_palette('magma', n_colors=(len(percs) + 2))
zorder = [11, 12, 13, 14, 15, 16, 17]
z = {p:z for p, z in zip(percs, zorder)}
c = {p:c for p, c in zip(percs, perc_cmap)}
grouped = kaki_data.groupby(['IPTGuM'])
df = pd.DataFrame([])
for g, d in grouped:
    for p in percs:     
        remainder = 100 - p
        low = remainder / 2
        upper = p + remainder / 2 
        _percs = np.percentile(d['fc_mu'], [low, upper])
        df = df.append({'percentile': p,
                         'IPTGuM': g,
                         'fc_low':_percs[0],
                         'fc_high': _percs[1]},
                      ignore_index=True)

for g, d in df.groupby('percentile'):
    ax7.fill_between(d['IPTGuM'], d['fc_low'], d['fc_high'], color=c[g], 
                    zorder=z[g], label=int(g))
  
# ##############################################################################
# EXPERIMENTAL MEAUREMENTS
# ##############################################################################
ax7.plot(data['IPTGuM'], data['fold_change'], '.', color='k', 
         markerfacecolor='w', label='data', zorder=1000)
leg = ax7.legend(title='percentile', loc='upper left', fontsize=6, ncol=2)
leg.get_title().set_fontsize(6)
plt.subplots_adjust(hspace=0.2, wspace=0.2)
plt.savefig('../../figures/FigSX_KaKi_only_post_predictive.pdf', 
            bbox_inches='tight')
