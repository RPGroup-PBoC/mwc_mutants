# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.viz
import mut.stats
import mut.thermo
mut.viz.plotting_style()
colors = mut.viz.color_selector('mut')
pboc_colors = mut.viz.color_selector('pboc')

# Load the data and the chains. 
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[(data['class'] == 'WT') | (data['class'] == 'IND')]
data = data[(data['mutant'] != 'Q294K') & (data['mutant'] != 'Q294R')]

# Load the chains from the MCMC.
chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_IND_strict.csv')
ind = np.argmax(chains['lnprobability'])

# Define the WT constants.
wt_ep = -13.9

# %%
# Set up the figure canvas.
c_range = np.logspace(-8, -2, 500)
fig, ax = plt.subplots(2, 2, figsize=(4.5, 4))
ax = ax.ravel()

# Format the axes.
for a in ax:
    a.set_xscale('log')
    a.tick_params(labelsize=8)

ax[1].set_xlim([1E-8, 1E-2])
ax[0].axis('off')
ax[1].set_xlabel('IPTG (M)', fontsize=8)
ax[1].set_ylabel('fold-change', fontsize=8)
ax[2].set_ylabel('dynamic range', fontsize=8)
ax[2].set_xlabel('$K_A / K_I$', fontsize=8)
ax[3].set_xlabel('$K_A$ (M)', fontsize=8)
ax[3].set_ylabel('$K_I$ (M)', fontsize=8)

# Add panel labels.
fig.text(-0.02, 0.95, '(A)', fontsize=8)
fig.text(0.48, 0.95, '(B)', fontsize=8)
fig.text(-0.02, 0.5, '(C)', fontsize=8)
fig.text(0.48, 0.5, '(D)', fontsize=8)

# Identify mutants and plot the true values. 
muts = data['mutant'].unique()
for i, m in enumerate(muts):

    # Get the most likely parameter estimates.
    ka_chain = chains['{}_Ka'.format(m)].values * 1E-6
    ki_chain = chains['{}_Ki'.format(m)].values * 1E-6
    ka_mode = ka_chain[ind]
    ki_mode = ki_chain[ind]

    # Compute the best-fit curve.
    fc_theo = mut.thermo.SimpleRepression(R=260, ep_r=wt_ep, ka=ka_mode, ki=ki_mode, 
                                        ep_ai=4.5, effector_conc=c_range).fold_change()

    cred_region = np.zeros((2, len(c_range)))
    for j, c in enumerate(c_range):
        arch = mut.thermo.SimpleRepression(R=260, ep_r=wt_ep, ka=ka_chain, ki=ki_chain, 
                                           ep_ai=4.5, effector_conc=c).fold_change()
        cred_region[:, j] = mut.stats.compute_hpd(arch, 0.95)

    # Plot the curves.    
    _ = ax[1].plot(c_range, fc_theo, lw=1.5, color=colors[m.upper()], label='__nolegend__')
    _ = ax[1].fill_between(c_range, cred_region[0, :], cred_region[1, :], color=colors[m.upper()], alpha=0.4)

# Plot the theoretical curves for the saturation ratio.     
def leak_kaki(R, ep_r, ep_ai=4.5, n_ns=4.6E6):
    return (1 + (1 + np.exp(-ep_ai))**-1 * (R / n_ns) * np.exp(-ep_r))**-1

def dyn_rng_kaki(R, ep_r, kaki, ep_ai=4.5, n=2, n_ns=4.6E6):
    leak = leak_kaki(R, ep_r, ep_ai, n_ns)
    pact = (1 + np.exp(-ep_ai) * (1 + kaki)**n)**-1
    return (1 + pact * (R / n_ns) * np.exp(-ep_r))**-1 - leak
kaki_range = np.logspace(1, 3, 500)
theo_dyn_rng = dyn_rng_kaki(260, wt_ep, kaki_range)
ax[3].set_yscale('log')
_ = ax[2].plot(kaki_range, theo_dyn_rng, 'k-', lw=1.5)

palette = sns.color_palette('Greys', n_colors = 6)
# Plot a contour of Ka Ki
ka_range = np.logspace(-5, -2, 500)
ki_range = np.logspace(-7, -5, 500)
KA, KI = np.meshgrid(ka_range, ki_range)
dyn_rng = mut.thermo.SimpleRepression(R=260, ep_r=wt_ep, ka=KA, ki=KI, effector_conc=1E9, ep_ai=4.5).dynamic_range()
_ = ax[3].contourf(ka_range, ki_range, dyn_rng, cmap='Blues_r', levels=[0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.98, 1])
contours = ax[3].contour(ka_range, ki_range, dyn_rng, colors='w', linewidths=0.75, levels=[0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.98, 1])
label_pos = [(1E-4, 8E-6), (1E-4, 2E-6), (5E-5, 8E-7), (3.7E-5, 2.3E-7), (5E-4, 1E-6), (5E-4, 5E-7)]
_ = ax[3].clabel(contours, fontsize=8, colors=palette, fmt='%0.2f', manual=label_pos)
ax[3].set_title('dynamic range', fontsize=8, backgroundcolor=pboc_colors['pale_yellow'], y=1.02)
# Plot the data.
grouped = data.groupby(['mutant', 'IPTGuM'])
for g, d in grouped:
    if g[0] == 'wt':
        face = 'w'
    else:
        face = colors[g[0]] 
    edge = colors[g[0].upper()]

    if g[1] == 0.0:
        label = g[0].upper()
    else:
        label = '__nolegend__'
    _ = ax[1].errorbar(g[1] / 1E6, d['fold_change'].mean(), d['fold_change'].std() / np.sqrt(len(d)), 
                       lw=1, color=edge, markerfacecolor=face, markeredgecolor=edge, markeredgewidth=1.5,
                       markersize=6, label=label, fmt='.')

# Plot the measured dynamic range.
grouped = data.groupby(['mutant'])
for g, d in grouped:
    ka = chains['{}_Ka'.format(g)] * 1E-6
    ki = chains['{}_Ki'.format(g)] * 1E-6
    kaki = ka / ki
    kaki_mode = kaki.values[ind]
    kaki_hpd = mut.stats.compute_hpd(kaki, 0.95)
    ka_hpd = mut.stats.compute_hpd(ka, 0.95)
    ki_hpd = mut.stats.compute_hpd(ki, 0.95)
    minlen = np.min([len(leak), len(sat)])
    leak = d[d['IPTGuM'] == 0]['fold_change']
    sat = d[d['IPTGuM'] == 5000]['fold_change']
    minlen = np.min([len(leak), len(sat)])
    dyn_rng = sat[:minlen-1].values - leak[:minlen-1].values
    

    if g == 'wt':
        face = 'w'
    else:
        face = colors[g]
    edge = colors[g.upper()]
    _ = ax[2].errorbar(kaki_mode, dyn_rng.mean(), np.std(dyn_rng)/np.sqrt(minlen), color=edge, 
    markerfacecolor=face, markeredgecolor=edge, lw=1.5, fmt='.', markersize=6, label=g.upper())
    _ = ax[2].hlines(dyn_rng.mean(), kaki_hpd[0], kaki_hpd[1], color=edge, lw=1, label='__nolegend__')
    _ = ax[3].plot(ka[ind], ki[ind], color=edge, markerfacecolor=face, markeredgecolor=edge, lw=1.5, marker='.', 
                    markersize=6, label=g.upper())
    _ = ax[3].hlines(ki[ind], ka_hpd[0], ka_hpd[1], color=edge, lw=1, label='__nolegend__')
    _ = ax[3].vlines(ka[ind], ki_hpd[0], ki_hpd[1], color=edge, lw=1, label='__nolegend__')

_ = ax[1].legend(loc='upper left', fontsize=8, handletextpad=0.3)
plt.tight_layout()
plt.savefig('../../figures/fig3_inducers.svg', bbox_inches='tight')
# %%
