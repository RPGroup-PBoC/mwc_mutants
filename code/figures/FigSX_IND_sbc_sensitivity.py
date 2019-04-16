# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
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
for i in range(2):
    ax[0, i].set_ylim([-5, 5])
    ax[0, i].set_xlim([-0.1, 1.1])
# ##############################################################################
# SHRINKAGE AND Z-SCORE
# ##############################################################################
for g, d in sbc_data.groupby(['model']):
    _ax = ax[0, axes[g]]
    ka = d[d['param']=='Ka']
    ki = d[d['param']=='Ki']
    if g != 'KaKi_only':
        ep = d[d['param']=='ep_AI']
        _ax.plot(ep['shrinkage'], ep['z_score'], ',', color='k')
    
    _ax.plot(ka['shrinkage'], ka['z_score'], ',', color=colors['red'])
    _ax.plot(ki['shrinkage'], ki['z_score'], ',', color=colors['blue'])
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
        _ax.step(ep_x, y, 'k')
    _ax.step(ka_x, y, color=colors['red']) 
    _ax.step(ki_x, y, color=colors['blue'])

   