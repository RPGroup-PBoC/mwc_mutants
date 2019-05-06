#-*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.viz
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()
sbc_data = pd.read_csv('../../data/csv/IND_sbc.csv')


# ##############################################################################
# Figure instantiation
# ##############################################################################
fig, ax = plt.subplots(1, 2, figsize=(7, 3))
for a in ax:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)
    a.set_xlim([0, 1.1])

# ##############################################################################
# KA_KI SHRNINKAGE PATHOLOGY
# ##############################################################################
for g, d in sbc_data.groupby(['model', 'sim_idx']):
    ka = d[d['param']=='Ka']
    ki = d[d['param']=='Ki']
    if ka['post_median'].values[0] >= ki['post_median'].values[0]:
        color=colors['red']
    else:
        color = 'k'
    if g[0]=='KaKi_only':
        ax[0].semilogy(ka['shrinkage'], ki['post_median'], '.', color=color)

# epai_data = sbc_data[sbc_data['param']=='ep_AI']
# ax[1].plot(epai_data['shrinkage'], epai_data['z'], 'k.')


d
