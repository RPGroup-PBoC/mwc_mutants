# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot 
import matplotlib.gridspec as gridspec
import mut.viz
import mut.thermo
import mut.stats
import scipy.stats
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()
ppc_data = pd.read_csv('../../data/csv/IND_prior_predictive_checks.csv')

ep_a_unique = ppc_data['ep_a'].unique()
ep_i_unique = ppc_data['ep_i'].unique()
ep_ai_unique = ppc_data[ppc_data['model']=='KaKi_epAI']['ep_ai'].unique()
# ##############################################################################
# FIGURE INSTANTIATION #
# ##############################################################################

fig = plt.figure(figsize=(6, 4))
gs = gridspec.GridSpec(4, 6)
ax1 = fig.add_subplot(gs[0:2, 0:2])
ax2 = fig.add_subplot(gs[0:2, 2:4])
ax3 = fig.add_subplot(gs[0:2, 4:])
ax4 = fig.add_subplot(gs[2:, 0:3])
ax5 = fig.add_subplot(gs[2:, 3:])
ax = [ax1, ax2, ax3, ax4, ax5]
for a in ax:
   a.xaxis.set_tick_params(labelsize=8) 
   a.yaxis.set_tick_params(labelsize=8)

# ##############################################################################
# TRUE PRIOR DISTRIBUTIONS
# ##############################################################################
ax1.plot(ep_a_unique, ep_i_unique, '.', ms=2, color=colors['red'])
ax2.plot(ep_ai_unique, ep_a_unique, '.', ms=2, color=colors['red'])
ax3.plot(ep_ai_unique, ep_i_unique, '.', ms=2, color=colors['red'])
# ax1.fill_between(epk_range, np.zeros(len(epa_epi_pdf)), epa_epi_pdf, color=colors['light_red'])
# ax2.plot(np.exp(epk_range), np.exp(epa_epi_pdf), '-', color=colors['red'])


