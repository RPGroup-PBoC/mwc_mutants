# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
sys.path.insert(0, '../../')
import mut.viz
import mut.stats
import mut.thermo
mut.viz.plotting_style()
colors = mut.viz.color_selector('mut')

# Load the compiled data.
data = pd.read_csv('../../data/csv/compiled_data.csv', comment='#')
data = data[(data['class'] == 'DNA') | (data['class'] == 'WT')]
data = data[data['fold_change'] > 0]


# Load the fitting chains.
chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_DNA_strict.csv')
ind = np.argmax(chains['lnprobability'].values)
modes = pd.DataFrame(chains.drop(
    columns=['lnprobability', 'sigma', 'Unnamed: 0']).iloc[ind]).reset_index()
modes.columns = ['parameter', 'mode']

# Define the architectural parameters.
Nns = 4.6E6
ka = 139E-6  # in M
ki = 0.53E-6  # in M
ep_r = {i.split('_')[0]: j for i, j in zip(modes['parameter'], modes['mode'])}


# Compute the leakiness, saturation, and dynamic range for each mutant.
grouped = data.groupby(['mutant', 'repressors'])
prop_df = pd.DataFrame([], columns=['mutant', 'repressors',  'leak_mean',
                                    'leak_sem', 'sat_mean', 'sat_sem', 'dyn_rng_mean', 'dyn_rng_sem', 'rkdna'])
for g, d in grouped:
    # Split by appropriate IPTG
    leak = d[d['IPTGuM'] == 0]['fold_change'].values
    sat = d[d['IPTGuM'] == 5000]['fold_change'].values
    min_length = np.min([len(leak), len(sat)])
    dyn_rng = sat[:min_length-1] - leak[:min_length-1]

    # Compute the mean and standard error for each one.
    leak_mean = leak.mean()
    leak_sem = np.std(leak) / np.sqrt(len(d))
    sat_mean = sat.mean()
    sat_sem = np.std(sat) / np.sqrt(len(d))
    dyn_rng_mean = dyn_rng.mean()
    dyn_rng_sem = np.std(dyn_rng) / np.sqrt(len(d))

    # Compute kdna/R
    rkda = (g[1] / Nns) * np.exp(-ep_r[g[0]])

    # Assemble the dictionary and append.
    prop_dict = {'mutant': g[0], 'repressors': g[1], 'leak_mean': leak_mean,
                 'leak_sem': leak_sem, 'sat_mean': sat_mean, 'sat_sem': sat_sem,
                 'dyn_rng_mean': dyn_rng_mean, 'dyn_rng_sem': dyn_rng_sem, 'rkdna': rkda}
    prop_df = prop_df.append(prop_dict, ignore_index=True)

# %%
# Define a function to compute the fold-change using rkdn
def fc_rkdna(rkdna, pact):
    return (1 + pact * rkdna)**-1


# Set up the MWC architecture and compute the extrema
mwc = mut.thermo.MWC(effector_conc=0, ka=ka, ki=ki, ep_ai=4.5)
leak_pact = mwc.leakiness()
sat_pact = mwc.saturation()

# Compute the theoretical curves.
rkdna_range = np.logspace(-3, 6, 500)
leak_theo = fc_rkdna(rkdna_range, leak_pact)
sat_theo = fc_rkdna(rkdna_range, sat_pact)
dyn_theo = sat_theo - leak_theo

# Set up the plots.
fig, ax = plt.subplots(2, 2, figsize=(5, 4))
ax = ax.ravel()
ax[-1].axis('off')
# Add panels.
fig.text(0, 1, '(A)', fontsize=8)
fig.text(0.5, 1, '(B)', fontsize=8)
fig.text(0, 0.5, '(C)', fontsize=8)

# Format the axes
for a in ax:
    a.set_xlabel(
        r'$\frac{R}{N_{ns}}e^{-\beta\Delta\varepsilon_{RA}}$', fontsize=8)
    a.set_xscale('log')
    a.set_xlim([1E-3, 1E6])
    a.set_ylim([-0.05, 1.1])
    a.tick_params(labelsize=8)

# Set the yaxes.
_ = ax[0].set_ylabel('leakiness', fontsize=8)
_ = ax[1].set_ylabel('saturation', fontsize=8)
_ = ax[2].set_ylabel('dynamic range', fontsize=8)
# Plot the theoretical curves in black.
_ = ax[0].plot(rkdna_range, leak_theo, 'k-', lw=1)
_ = ax[1].plot(rkdna_range, sat_theo, 'k-', lw=1)
_ = ax[2].plot(rkdna_range, dyn_theo, 'k-', lw=1, label='prediction')

# Plot the data.
grouped = prop_df.groupby(['mutant'])
for g, d in grouped:
    if g == 'wt':
        face = 'w'
        edge = 'k'
    else:
        face = colors[g]
        edge = colors[g]
    _ = ax[0].errorbar(d['rkdna'], d['leak_mean'], d['leak_sem'], fmt='.', ms=5, color=colors[g.upper()],
                       lw=0.75, markerfacecolor=face, markeredgecolor=edge, markeredgewidth=0.75,
                       label=g)
    _ = ax[1].errorbar(d['rkdna'], d['sat_mean'], d['sat_sem'], fmt='.', ms=5, color=colors[g.upper()],
                       lw=0.75, markerfacecolor=face, markeredgecolor=edge, markeredgewidth=0.75,
                       label=g)
    _ = ax[2].errorbar(d['rkdna'], d['dyn_rng_mean'], d['dyn_rng_sem'], fmt='.', color=colors[g.upper()],
                       ms=5, lw=0.75, markerfacecolor=face, markeredgecolor=edge, markeredgewidth=0.75,
                       label=g)
ax[2].legend(bbox_to_anchor=(2, 1.0), fontsize=8)

plt.tight_layout()
plt.savefig('../../figures/fig5.pdf', bbox_inches='tight')
