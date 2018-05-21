# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.viz
import mut.stats
import mut.thermo
mut.viz.plotting_style()
colors = mut.viz.color_selector('mut')
pboc_colors = mut.viz.color_selector('pboc')

# Load the data and chains.
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[(data['class'] == 'DNA') | (data['class'] == 'WT')]
chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_DNA_strict.csv')

# Define the constants.
wt_ka = 139E-6  # M
wt_ki = 0.53E-6  # M
muts = data['mutant'].unique()
ind = np.argmax(chains['lnprobability'])
c_range = np.logspace(-8, -2, 500)

# Set up a function to compute the saturation and leakiness with rkdna.
# %%


def leak_rkdna(rkdna, ep_ai=4.5, n_ns=4.6E6):
    pact = (1 + np.exp(-ep_ai))**-1
    return (1 + pact * rkdna)**-1


def sat_rkdna(rkdna, kaki, n=2, ep_ai=4.5, n_ns=4.6E6):
    pact = (1 + np.exp(-ep_ai) * kaki**n)**-1
    return (1 + pact * rkdna)**-1


# Compute the theoretical curves
rkdna_range = np.logspace(-3, 5, 500)
leak = leak_rkdna(rkdna_range)
sat = sat_rkdna(rkdna_range, wt_ka/wt_ki)
dyn_rng = sat - leak

# Set up the figure canvas.
fig, ax = plt.subplots(2, 2, figsize=(4.5, 4))
ax = ax.ravel()
for a in ax:
    a.set_xscale('log')
    a.tick_params(labelsize=8)
    a.set_ylim([-0.15, 1.15])
# Add appropriate labels. ax[1].set_ylabel('fold-change', fontsize=8) ax[2].set_ylabel('leakiness', fontsize=8)
ax[0].axis('off')
ax[1].set_ylabel('fold-change', fontsize=8)
ax[3].set_ylabel('dynamic range', fontsize=8)
ax[2].set_ylabel('leakiness', fontsize=8)
ax[1].set_xlabel('IPTG (M)', fontsize=8)
ax[2].set_xlabel(r'$\frac{R}{N_{NS}}e^{-\beta\Delta\varepsilon_{RA}}$')
ax[3].set_xlabel(r'$\frac{R}{N_{NS}}e^{-\beta\Delta\varepsilon_{RA}}$')
ax[1].set_xlim([1E-8, 1E-2])
# Add panel labels.
fig.text(0, 0.95, '(A)', fontsize=8)
fig.text(0.5, 0.95, '(B)', fontsize=8)
fig.text(0, 0.5, '(C)', fontsize=8)
fig.text(0.5, 0.5, '(D)', fontsize=8)

# Plot the theoretical curves for the rkdna plots.
_ = ax[2].plot(rkdna_range, leak, 'k', lw=1.5)
_ = ax[3].plot(rkdna_range, dyn_rng, 'k', lw=1.5)
for a in ax[2:]:
    a.set_xlim([1E-3, 1E5])

# Plot the theotretical titration curves.
for i, m in enumerate(muts):

    # Get the sampling chains.
    epr_mode = chains.iloc[ind]['{}_eps_r'.format(m)]
    epr_chains = chains['{}_eps_r'.format(m)]

    # Compute the theoretical fold-change.
    fc_theo = mut.thermo.SimpleRepression(R=260, ep_r=epr_mode, ka=wt_ka, ki=wt_ki, ep_ai=4.5,
                                          effector_conc=c_range).fold_change()

    # Compute the credible region.
    cred_region = np.zeros((2, len(c_range)))
    for j, c in enumerate(c_range):
        arch = mut.thermo.SimpleRepression(R=260, ep_r=epr_chains, ka=wt_ka, ki=wt_ki, ep_ai=4.5,
                                           effector_conc=c).fold_change()
        cred_region[:, j] = mut.stats.compute_hpd(arch, 0.95)

    # Plot the theoretical curves.
    _ = ax[1].plot(c_range, fc_theo, color=colors[m.upper()],
                   lw=1.5, label='__nolegend__')
    _ = ax[1].fill_between(c_range, cred_region[0, :], cred_region[1, :], color=colors[m.upper()], alpha=0.5)


# Plot the titration data.
grouped = data[data['repressors'] == 260].groupby(['mutant', 'IPTGuM'])
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
                       color=edge, lw=1, markeredgecolor=edge, markerfacecolor=face, markeredgewidth=1, markersize=6, label=label,
                       fmt='.')


# Plot the saturation and dynamic range.
markers = ['o', 's', 'D', '^']
reps = np.sort(data['repressors'].unique())
marker_dict = {i: j for i, j in zip(reps, markers)}
grouped = data.groupby(['mutant', 'repressors'])
for g, d in grouped:
    if g[0] == 'wt':
        face = 'w'
    else:
        face = colors[g[0]]
    edge = colors[g[0].upper()]
    marker = marker_dict[g[1]]

    # Compute rkDNA.
    epr_mode = chains.iloc[ind]['{}_eps_r'.format(g[0])]
    epr_hpd = mut.stats.compute_hpd(chains['{}_eps_r'.format(g[0])], 0.95)

    rkdna_mode = (g[1]/4.6E6) * np.exp(-epr_mode)
    rkdna_min = (g[1]/4.6E6) * np.exp(-epr_hpd[0])
    rkdna_max = (g[1]/4.6E6) * np.exp(-epr_hpd[1])

    # Find the leakiness and saturation
    leak = d[d['IPTGuM'] == 0]['fold_change'].values
    sat = d[d['IPTGuM'] == 5000]['fold_change'].values
    min_len = np.min([len(leak), len(sat)])
    dyn_rng = sat[:min_len-1] - leak[:min_len-1]

    _ = ax[2].errorbar(rkdna_mode, leak.mean(), yerr=leak.std() / np.sqrt(min_len), color=edge,
                   markerfacecolor=face, markeredgecolor=edge, markeredgewidth=1.5, markersize=4, fmt=marker,
                   label='__nolegend__')
    _ = ax[2].hlines(leak.mean(), rkdna_min, rkdna_max, color=edge, lw=1, label='__nolegend__')


    _ = ax[3].errorbar(rkdna_mode, dyn_rng.mean(), yerr=dyn_rng.std() / np.sqrt(min_len), color=edge,
                   markerfacecolor=face, markeredgecolor=edge, markeredgewidth=1.5, markersize=4, fmt=marker,
                   label='__nolegend__')
    _ = ax[3].hlines(dyn_rng.mean(), rkdna_min, rkdna_max, color=edge, lw=1, label='__nolegend__')

# Add the apprpriate legends.
for m in marker_dict:
    _ = ax[2].plot([], [], color='slategray', linestyle='none', marker=marker_dict[m], markersize=4, alpha=0.5, label=int(m))
leg = ax[2].legend(loc='upper right', fontsize=8, title='rep. / cell', handletextpad=0.01)

ax[1].legend(loc='upper left', handletextpad=0.005, fontsize=8, labelspacing=0.07, borderpad=0.02)
leg.get_title().set_fontsize(8)
plt.tight_layout()
plt.savefig('../../figures/fig2_DNA.svg', bbox_inches='tight')
# %%
dyn_rng