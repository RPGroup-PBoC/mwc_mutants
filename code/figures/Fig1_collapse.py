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
constants = mut.thermo.load_constants()
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()
colors = {c:i for c, i in colors.items() if ('light' not in c) & ('pale' not in c)}

# Load the data.
old_gods = pd.read_csv('../../data/csv/Garcia2011_Brewster2014_data.csv')
new_gods = pd.read_csv('../../data/csv/RazoMejia2018_data.csv')
new_gods = new_gods[new_gods['repressors'] > 0]

# Prune non-physical fold-change. 
old_gods = old_gods[(old_gods['fold_change'] >= -0.1) & (old_gods['fold_change'] <= 1.2)]
# Summarize the induction data. 
new_gods.loc[:, 'fold_change'] = new_gods['fold_change_A']
new_gods = pd.DataFrame(new_gods.groupby(
    ['repressors', 'operator', 'IPTG_uM']).apply(mut.stats.compute_mean_sem)).reset_index()


# Plot the collapse curve. 
bohr_range = np.linspace(-10, 10, 200)
collapse = (1 + np.exp(-bohr_range))**-1

# Instantiate the figure.
fig, ax = plt.subplots(1, 1, figsize=(3, 2))
ax.xaxis.set_tick_params(labelsize=8)
ax.yaxis.set_tick_params(labelsize=8)
ax.set_xlabel('Bohr parameter [$k_BT$]', fontsize=8)
ax.set_ylabel('fold-change', fontsize=8)

# Plot the master curve. 
_ = ax.plot(bohr_range, collapse, 'k-', label='master curve')

# Iterate through the brewster and garcia data and compute the bohr parameter.
authors = {'garcia':colors['red'], 'brewster':colors['purple']}

for g, d in old_gods.groupby(['author', 'energy', 'repressor']):
    # Compute and plot bohrconstants.
    bohr = mut.thermo.SimpleRepression(R=g[2], ep_r=g[1], ka=constants['Ka'],
                                      ki=constants['Ki'], effector_conc=0, n_sites=constants['n_sites'], 
                                       n_ns=constants['Nns'], ep_ai=constants['ep_AI']).bohr_parameter()
    _ = ax.plot(bohr, d['fold_change'], 'o', color=authors[g[0]], ms=1.5, label='__nolegend__', alpha=0.75)

# Plot the flow data. 
for g, d in new_gods.groupby(['operator', 'IPTG_uM', 'repressors']):
    bohr = mut.thermo.SimpleRepression(R=2 * g[2], ep_r=constants[g[0]], ka=constants['Ka'],
                                      ki=constants['Ki'], ep_ai=constants['ep_AI'],
                                      n_sites=constants['n_sites'], n_ns=constants['Nns'],
                                      effector_conc=g[1]).bohr_parameter()
    _ = ax.plot(bohr, d['mean'], 'o', color=colors['blue'], ms=1.5, label='__nolegend__', alpha=0.75)
    
# Add the appropriate labels
ax.plot([], [], 'o', ms=2, color=colors['red'], label='Garcia & Phillips\n2011')
ax.plot([], [], 'o', ms=2, color=colors['purple'], label='Brewster et al.\n2014')
ax.plot([], [], 'o', ms=2, color=colors['blue'], label='Razo-Mejia et al.\n2018')

ax.legend(fontsize=6, handletextpad=0.1, handlelength=1, loc='upper left')
plt.tight_layout()
plt.savefig('Fig1_collapse.svg')

