#%%
# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.thermo
import mut.stats
import mut.viz
colors = mut.viz.pub_style()

# Define the experimental constants.
n_ns = 4.6E6
epR = -13.9
ep_ai = 4.5
R = 260

# Load the data and restrict to inducer mutants.
data = pd.read_csv('../../data/csv/compiled_data.csv')
IND = data[(data['class']=='IND') | (data['class'] == 'WT')]
muts = ['F164T', 'Q294V', 'Q294R']
# Load the fit statistics and the flat chains.
samples = pd.read_csv('../../data/mcmc/IND_O2_KaKi_fit_chains.csv')
stats = pd.read_csv('../../data/mcmc/IND_O2_KaKi_fit_statistics.csv')


#%% 
# Define the color key. 
color_key = {'F164T':colors['red'], 'Q294V':colors['green'], 'Q294R':colors['blue']}

#Instantiate the figure canvas
fig, ax = plt.subplots(1, 1,  figsize=(6,4))
ax.set_xscale('log')
ax.set_xlabel('IPTG [M]')
ax.set_ylabel('fold-change')

# Plot the fits
c_range = np.logspace(-8, -2, 500)
for i, m in enumerate(muts):
    # Extract the parameters and flatchains.
    ka = stats[stats['parameter']=='ka_{}'.format(m)]['mode'].values[0]
    ki = stats[stats['parameter']=='ki_{}'.format(m)]['mode'].values[0]
    ka_chain = samples['ka_{}'.format(m)]
    ki_chain = samples['ki_{}'.format(m)]

    # Comput the foldchange.
    arch = mut.thermo.SimpleRepression(R=260, ka=ka/1E6, ki=ki/1E6, ep_r=epR, effector_conc=c_range, ep_ai=4.5)
    fc_theo = arch.fold_change()
    _ = ax.plot(c_range, fc_theo, '-', lw=1, color=color_key[m])

_grouped = IND.groupby(['mutant', 'IPTGuM']).apply(mut.stats.compute_mean_sem)
grouped_df = pd.DataFrame(_grouped).reset_index()
grouped = grouped_df.groupby(['mutant', 'IPTGuM'])
for g, d in grouped:
    if g[0] != 'wt':    
        _ = ax.errorbar(g[1] / 1E6, d['mean'], d['sem'], fmt='o',
        color=color_key[g[0]])

stats[stats['parameter']=='ka_{}'.format(m)]['mode'].values[0]