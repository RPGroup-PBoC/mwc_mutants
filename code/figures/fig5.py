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
    leak_sem = np.std(leak) / np.sqrt(len(leak))
    sat_mean = sat.mean()
    sat_sem = np.std(sat) / np.sqrt(len(sat))
    dyn_rng_mean = dyn_rng.mean()
    dyn_rng_sem = np.std(dyn_rng) / np.sqrt(len(leak))

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

# Instantiate the figure
c_range = np.logspace(-8, -2, 300)
fig, ax  = plt.subplots(2, 2, figsize=(6, 4.5))
fig.text(0, 0.95, '(A)', fontsize=8)
fig.text(0.5, 0.95, '(B)', fontsize=8)
fig.text(0, 0.5, '(C)', fontsize=8)
fig.text(0.5, 0.5, '(D)', fontsize=8)
ax = ax.ravel()
ax[0].axis('off')
for a in ax:
    a.tick_params(labelsize=8)
    a.set_xscale('log')
    a.set_ylim([-0.12, 1.12])

# Format the axes.
ax[1].set_xlabel('IPTG (M)', fontsize=8)
ax[1].set_ylabel('fold-change', fontsize=8)
ax[2].set_xlabel(r'$\frac{R}{N_{NS}}e^{-\beta\Delta\varepsilon_{RA}}$', fontsize=8)
ax[2].set_ylabel('leakiness', fontsize=8)
ax[3].set_xlabel(r'$\frac{R}{N_{NS}}e^{-\beta\Delta\varepsilon_{RA}}$', fontsize=8)
ax[3].set_ylabel('dynamic range', fontsize=8)
ax[1].set_xlim([1E-8, 1E-2])
for a in ax[2:]:
    a.set_xlim([1E-3, 1E6])

# Plot the fits and the credible regions
mutants = data['mutant'].unique()
for i, m in enumerate(mutants):
 # Get the parameter values
    ep_r = modes[modes['parameter']=='{}_eps_r'.format(m)]['mode'].values[0]
    
    fc = mut.thermo.SimpleRepression(R=260, ep_r=ep_r, ka=ka, ki=ki, ep_ai=4.5,
                                            effector_conc=c_range).fold_change()
    cred_region = np.zeros((2, len(c_range)))
    for j, c in enumerate(c_range):
        prob = mut.thermo.SimpleRepression(R=260, ep_r=chains['{}_eps_r'.format(m)], ka=ka, ki=ki, ep_ai=4.5,
                                            effector_conc=c).fold_change()
        cred_region[:, j] = mut.stats.compute_hpd(prob, 0.95) 

    # Plot the fits and credible regions.
    _ = ax[1].plot(c_range, fc, color=colors[m.upper()], lw=1.5)
    _ = ax[1].fill_between(c_range, cred_region[0, : ], cred_region[1, :], color=colors[m.upper()],
                            alpha=0.5)


# Plot the fold-change data
_data = data[data['repressors'] == 260].copy()
grouped = _data.groupby(['mutant', 'IPTGuM'])
for g, d in grouped:
    if g[0] == 'wt':
        face = 'w'
    else:
        face = colors[g[0].upper()]
   
    if g[1] == 0:
        if g[0] == 'wt':
            label = 'wild-type'
        else:
            label = g[0].upper()
    else:
        label = '__nolegend__'
    ax[1].errorbar(g[1]/1E6, d['fold_change'].mean(), d['fold_change'].std() / np.sqrt(len(d)), 
        color=colors[g[0].upper()], markeredgecolor=colors[g[0].upper()], markerfacecolor=face, 
        markeredgewidth=1.5, fmt='.', ms=6, lw=1, linestyle='none', label=label)


markers = ['o', 's', 'D', '^']
reps = prop_df['repressors'].unique() 
marker_dict = {i:j for i, j in zip(reps, markers)}
# Plot the leakiness and saturation.
grouped = prop_df.groupby(['mutant', 'repressors'])
for g, d in grouped:
    if g[0] == 'wt':
        face = 'w'
    else:
        face = colors[g[0].upper()]
    ax[2].errorbar(d['rkdna'], d['leak_mean'], d['leak_sem'], color=colors[g[0].upper()], fmt=marker_dict[g[1]],
    markerfacecolor=face, markeredgecolor=colors[g[0].upper()], ms=4, markeredgewidth=1.5, label='__nolegend__')
    ax[3].errorbar(d['rkdna'], d['dyn_rng_mean'], d['dyn_rng_sem'], color=colors[g[0].upper()], fmt=marker_dict[g[1]],
    markerfacecolor=face, markeredgecolor=colors[g[0].upper()], ms=4, markeredgewidth=1.5)

for r, m in marker_dict.items():
    ax[2].plot([], [], color='slategray', ms=4, marker=m, label=int(r), linestyle='none', alpha=0.5)
ax[1].legend(fontsize=8, loc='upper left', labelspacing=0.3, handletextpad=0.2)
ax[1].text(0.8, 0.15, '$R$ = 260', fontsize=8, transform=ax[1].transAxes)
leg = ax[2].legend(title='rep. / cell', fontsize=8)
leg.get_title().set_fontsize(8)

# Plot the theoretical values for leakiness and dynamic range
ax[2].plot(rkdna_range, leak_theo, 'k-', lw=1.5)
ax[3].plot(rkdna_range, dyn_theo, 'k-', lw=1.5)
plt.tight_layout()
plt.savefig('../../figures/fig3_plots.svg')








# %% Set up the plots.
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
