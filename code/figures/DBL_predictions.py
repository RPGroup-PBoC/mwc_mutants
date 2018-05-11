#%%
# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sys.path.insert(0, '../../')
import mut.thermo
import mut.stats
import mut.viz
colors = mut.viz.pub_style()

# Load the statistics from each mutant class
DNA_stats = pd.read_csv('../../data/mcmc/DNA_O2_epR_fit_statistics.csv')
IND_stats = pd.read_csv('../../data/mcmc/IND_O2_KaKi_fit_statistics.csv')

# Load the data and restrict to the double mutants.
data = pd.read_csv('../../data/csv/compiled_data.csv')
DBL = data[data['class']=='DBL']
DBL = DBL[DBL['fold_change'] >= 0]

# Get the  different mutant IDS.
DNA_muts = ['Y20I', 'Q21M', 'Q21A']
IND_muts = ['F164T', 'Q294V']

# Set up the figure canvas.
fig, ax = plt.subplots(1, 1, figsize=(6, 4))

ax.set_xscale('log')

_grouped = DBL.groupby(['mutant', 'IPTGuM']).apply(mut.stats.compute_mean_sem)
_df = pd.DataFrame(_grouped).reset_index()
grouped = _df.groupby(['mutant'])

# Iterate through each mutant and plot the data.
color_key = {'F164T': colors['red'], 'Q294V':colors['blue']}
# axes = {'Q21M': ax[0], 'Q21A': ax[1], 'Y20I':ax[2]}
c_range = np.logspace(-8, -2, 500)
color_palette = sns.color_palette('deep')
i = 0
for g, d in grouped:
    if g.split('-')[1] != 'Q294K':

        # Plot the best fit
        _dna, _ind = g.split('-')

        epR = DNA_stats[DNA_stats['parameter']=='ep_R_{}'.format(_dna)]['mode'].values[0]
        ka = IND_stats[IND_stats['parameter']=='ka_{}'.format(_ind)]['mode'].values[0]
        ki = IND_stats[IND_stats['parameter']=='ki_{}'.format(_ind)]['mode'].values[0]
        arch = mut.thermo.SimpleRepression(ep_r=epR, R=260, effector_conc=c_range, ka=ka/1E6, ki=ki/1E6,
        ep_ai=4.5, n_sites=2)        
        fc_theo = arch.fold_change()
        _ = ax.errorbar(d['IPTGuM']/1E6, d['mean'], d['sem'], fmt='o', linewidth=1, label=g,
         ms=2, color=color_palette[i])
        _ = ax.plot(c_range, fc_theo, lw=1, color=color_palette[i])
        i += 1
ax.legend()
ax.set_xlabel('IPTG [M]', fontsize=10)
ax.set_ylabel('fold-change', fontsize=10)
plt.tight_layout()
plt.savefig('../../figures/DBL_predictions.pdf', bbox_inches='tight')




