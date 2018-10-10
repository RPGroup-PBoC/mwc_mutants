# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
sys.path.insert(0, '../../')
import mut.bayes
import mut.thermo
import mut.stats
import mut.viz
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()
bright_colors = {m:i  for m, i in colors.items() if '_' not in m}

# Load the summarized data and Hill fit samples.
data = pd.read_csv('../../data/csv/summarized_data.csv')
hill_samples = pd.read_csv('../../data/csv/FigS1_Hill_samples.csv')
hill_stats = pd.read_csv('../../data/csv/FigS1_Hill_stats.csv')

# Screen the data. 
data = data[(data['class']=='DNA') & (data['operator']=='O2')]

# Extract the modes from the hill function parameters. 
modes = hill_samples.iloc[np.argmax(hill_samples['logp'])]

# Define the colors
rep_colors = {r: list(bright_colors.values())[i] for i, r in enumerate(data['repressors'].unique())}

# Define a function to compute the hill. 
def generalized_hill(a, b, c, k, n):
    return a + b * ((c/k)**n / (1 + (c/k)**n))


# Instantiate the figure canvas. 
fig = plt.figure(figsize=(6, 6.5))
gs = gridspec.GridSpec(11, 3)

# Properly format each axis
axes = []
subaxes = []
for i in range(3):
    _ax = fig.add_subplot(gs[:3, i]) 
    _ax.set_xscale('log')
    axes.append(_ax)
    _ax1 = [fig.add_subplot(gs[4 + j: 4+j+2, i]) for j in range(0, 8, 2)]
    for j, a in enumerate(_ax1):
        if j == 3:
            a.set_ylim([0, 2])
            a.set_xticks([0, 1, 2, 3])
            a.set_xticklabels(data['repressors'].unique().astype(int), rotation='vertical')
        else: 
            a.set_xticklabels(['', '', '', '']) 
        a.set_xlim([-0.5, 3.5])
        a.xaxis.set_tick_params(labelsize=8)
        a.yaxis.set_tick_params(labelsize=8)
    
    subaxes.append(_ax1) 
axes[0].set_ylabel('fold-change', fontsize=8)
# Properly format each subaxis
mut_ax = {'Y20I':axes[0], 'Q21A':axes[1], 'Q21M':axes[2]}
for t, a in mut_ax.items():
    a.set_title(t, backgroundcolor=colors['pale_yellow'], y=1.04, fontsize=8)
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_xlabel('IPTG [M]', fontsize=8)
    a.set_ylim([0, 1])
    a.set_xlim([1E-8, 1E-2])
    
    if t != 'Y20I':
        a.set_yticklabels(['', '', '', ''])
    
mut_subax = {'Y20I':subaxes[0], 'Q21A':subaxes[1], 'Q21M':subaxes[2]}
c_range = np.logspace(-2, 4, 200)
rep_pos = {r:i for i, r in enumerate(data['repressors'].unique())}

pars = ['a', 'b', 'K', 'n']

# Group and plot the data. 
grouped = data.groupby(['mutant', 'repressors'])
for g, d in grouped:
    mut_ax[g[0]].errorbar(d['IPTGuM'] / 1E6, d['mean'], d['sem'], fmt='o', lw=1, linestyle='none', ms=3,
                         label=int(g[1]), color=rep_colors[g[1]])    
    # Compute the line of best fit. 
    idx = f'{g[0]}.{int(g[1])}'
    hill_mode = generalized_hill(modes[f'a.{idx}'], modes[f'b.{idx}'],
                                c_range, modes[f'K.{idx}'], modes[f'n.{idx}'])
    mut_ax[g[0]].plot(c_range / 1E6, hill_mode, lw=1, color=rep_colors[g[1]])
    for j, p in enumerate(pars):
        mut_subax[g[0]][j].plot(rep_pos[g[1]], modes[f'{p}.{idx}'], 'o', markerfacecolor='w', markersize=4, 
                             markeredgecolor=rep_colors[g[1]])
        if p == 'K': 
            label = '$K$ [µM]'
        else:
            label = '$%s$' %p
            
        mut_subax[g[0]][j].set_ylabel(label, fontsize=8)
    mut_subax[g[0]][-1].set_xlabel('repressors / cell', fontsize=8)

# Hardcode the bounds for  K
mut_subax['Y20I'][2].set_ylim([0,  25])
mut_subax['Q21A'][2].set_ylim([0,  50])
mut_subax['Q21M'][2].set_ylim([50, 190])

# Adjust the bounds for a and b
for m in data['mutant'].unique():
    for i in range(2):
        mut_subax[m][i].set_ylim([0,1])
    if m == 'Q21M':
        mut_subax[m][0].set_ylim([1E-3, 2E-2])
        
# Add the legend where necessary. 
_leg = mut_ax['Q21M'].legend(loc='upper left', fontsize=7, handletextpad=0.15, title='rep. / cell')
_leg.get_title().set_fontsize(7)
plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.savefig('hill_fits.svg')


