# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.viz
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load the SBC statistics
sbc_data = pd.read_csv('../../data/csv/IND_sbc.csv')

# Reform the dataframe to only focus on Ka and Ki
df = pd.DataFrame([])
for g, d in sbc_data.groupby(['sim_idx', 'model']):
    ka = d[d['param']=='Ka'][['post_median', 'shrinkage']].values[0]
    ki = d[d['param']=='Ki'][['post_median', 'shrinkage']].values[0]
    df = df.append({'Ka_shrinkage':ka[1], 'Ki_shrinkage':ki[1], 
                    'Ka_median':ka[0], 'Ki_median':ki[0], 
                    'model':g[1], 'sim':g[0], 'idx':ka[1] <= ki[1]},
                    ignore_index=True)

# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(2, 1, figsize=(6, 3))

# Define the axes
axes = {'KaKi_only':ax[0], 'KaKi_epAI':ax[1]}

# Adjust scaling and set labels
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_xlim([-1, 1.1])
    a.set_yscale('log')
    a.set_xlabel('shrinkage', fontsize=8)
    a.set_ylabel('$K_A$ [µM]', fontsize=8)
ax[0].set_title('$K_A$ and $K_I$ modified', fontsize=8)
ax[0].set_title(r'$K_A$, $K_I$, and $\Delta\varepsilon_{AI}$ modified', fontsize=8)

# ##############################################################################
# ECDFS
# ##############################################################################
for g, d in df.groupby('model'):
    gt_s = d[d['idx']==0]['Ka_shrinkage'].values
    gt_v = d[d['idx']==0]['Ka_median'].values
    lt_s = d[d['idx']==1]['Ka_shrinkage'].values
    lt_v = d[d['idx']==1]['Ka_median'].values
    pool_s = d[d['idx']==1]['Ka_shrinkage'].values
    pool_v = d[d['idx']==1]['Ka_median'].values
    _ax = axes[g]
    _ax.plot(gt_s, gt_v, '.', color=colors['red'], 
             label='$K_A > K_I$', alpha=0.5)
    _ax.plot(lt_s, lt_v, '.', color=colors['blue'], 
             label='$K_A < K_I$', alpha=0.5)

for a in ax:
    a.legend(fontsize=8, ncol=2, loc='lower left')
plt.tight_layout()