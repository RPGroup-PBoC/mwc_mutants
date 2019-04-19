# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot  as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import mut.viz
import mut.thermo
import mut.stats
import scipy.stats
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()
ppc_data = pd.read_csv('../../data/csv/IND_prior_predictive_checks.csv')

ep_a_unique = ppc_data['ka'].unique()
ep_i_unique = ppc_data['ki'].unique()
ep_ai_unique = ppc_data[ppc_data['model']=='KaKi_epAI']['ep_ai'].unique()

# ##############################################################################
# FIGURE INSTANTIATION #
# ##############################################################################
fig, ax = plt.subplots(2, 3, figsize=(6, 4))
ax1, ax2, ax3, ax4, _, ax5 = ax.ravel()
for a in ax.ravel(0):
   a.xaxis.set_tick_params(labelsize=6) 
   a.yaxis.set_tick_params(labelsize=6)
_.set_axis_off()

# Add axis labels
ax1.set_xlabel(r'$K_A$ [µM]', fontsize=8, labelpad=0.1)
ax1.set_ylabel(r'$K_I$ [µM]', fontsize=8, labelpad=0.1)
ax2.set_xlabel(r'$\Delta\varepsilon_{AI}$ [$k_BT$]', fontsize=8, labelpad=0.1)
ax2.set_ylabel(r'$K_A$ [µM]', fontsize=8, labelpad=0.1)
ax3.set_xlabel(r'$\Delta\varepsilon_{AI}$ [$k_BT$]', fontsize=8)
ax3.set_ylabel(r'$K_I$ [µM]', fontsize=8, labelpad=0.1)

# Set scaling
ax4.set_xscale('symlog', linthreshx=1E-2)
ax5.set_xscale('symlog', linthreshx=1E-2)

# Set limits
models = ['$K_A$ and $K_I$ only', r'$K_A$, $K_I$, and $\Delta\varepsilon_{AI}$']
for i, a in enumerate([ax4, ax5]):
   a.set_xlim([-0.001, 5E3])
   a.set_yticks([-0.25, 0, 0.25, 0.5, 0.75, 1, 1.25])
   a.set_xlabel('IPTG [µM]', fontsize=8)
   a.set_ylabel('fold-change', fontsize=8)
   a.set_title(models[i], fontsize=8, y=1.08, backgroundcolor=colors['pale_yellow'])
   a.set_xticks([0, 1E-2, 1E0, 1E2, 1E4])

# Add panel labels
fig.text(0, 0.90, '(A)', fontsize=8)
fig.text(0, 0.45, '(B)', fontsize=8)

# Define the axes
axes = {'KaKi_only':  ax4, 'KaKi_epAI': ax5}
# ##############################################################################
# SAMPLED PRIOR DISTRIBUTIONS
# ##############################################################################
ax1.loglog(ep_a_unique, ep_i_unique, '.', ms=2, color=colors['red'])
ax2.semilogy(ep_ai_unique, ep_a_unique, '.', ms=2, color=colors['blue'])
ax3.semilogy(ep_ai_unique, ep_i_unique, '.', ms=2, color=colors['blue'])


# ##############################################################################
#  PRIOR PREDICTIVE CHECKS
# ##############################################################################
percs = [99, 95, 80, 50, 20, 10, 5]
cmap_kakionly = {p:c for p, c in zip(percs, sns.color_palette('Reds', len(percs)))}
cmap_kakiepai = {p:c for p, c, in zip(percs, sns.color_palette('Blues', len(percs)))}
zorder = {p:i for p, i in zip(percs, [10, 11, 12, 13, 14, 15, 16, 17])}

# Compute the percentiles of the simulations. 
grouped = ppc_data.groupby(['IPTGuM', 'model'])
df = pd.DataFrame([], columns=['percentile', 'IPTGuM', 'fc_low', 'fc_high', 'model'])
for g, d in grouped:
    for p in percs:     
        remainder = 100 - p
        low = remainder / 2
        upper = p + remainder / 2 
        _percs = np.percentile(d['fc_draw'], [low, upper])
        df = df.append({'percentile': p,
                         'IPTGuM': g[0],
                         'fc_low':_percs[0],
                         'fc_high': _percs[1],
                         'model': g[1]},
                      ignore_index=True)

for g, d in  df.groupby(['model', 'percentile']):
   _ax = axes[g[0]]
   if g[0] == 'KaKi_only':
      cmap = cmap_kakionly
   else:
      cmap = cmap_kakiepai
   _ax.fill_between(d['IPTGuM'], d['fc_low'], d['fc_high'], color=cmap[g[1]],
                  zorder=zorder[g[1]], label = g[1])
leg = ax4.legend(title='percentile', fontsize=6, bbox_to_anchor=(1.9, 1))
leg.get_title().set_fontsize(6)
leg = ax5.legend(title='percentile', fontsize=6, bbox_to_anchor=(-0.7, 1))
leg.get_title().set_fontsize(6)
plt.subplots_adjust(wspace=0.6, hspace=0.6)
plt.savefig('../../figures/FigSX_IND_prior_predictive_checks.pdf', 
         bbox_inches='tight', background='white')


