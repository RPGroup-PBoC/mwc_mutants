# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
sys.path.insert(0, '../../')
import mut.viz
import mut.stats
import mut.thermo
mut.viz.plotting_style()
colors = mut.viz.color_selector('mut')
pboc_colors = mut.viz.color_selector('pboc')

# Load the data.
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[(data['class'] == 'DBL') | (data['class'] == 'WT')]
dbls = ['Y20I-F164T', 'Y20I-Q294V', 'Q21A-F164T',
        'Q21A-Q294V', 'Q21M-F164T', 'Q21M-Q294V']
dbls_labels = ['Y20I\nF164T', 'Y20I\nQ294V', 'Q21A\nF164T',
               'Q21A\nQ294V', 'Q21M\nF164T', 'Q21M\nQ294V']

# %% Load the chains.
DNA_chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_DNA_strict.csv')
IND_chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_IND_strict.csv')
DBL_chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_DBL_strict.csv')
DNA_ind = np.argmax(DNA_chains['lnprobability'].values)
IND_ind = np.argmax(IND_chains['lnprobability'].values)
DBL_ind = np.argmax(DBL_chains['lnprobability'].values)

# %% Compute the theoretical collapse.
bohr_param = np.linspace(-10, 10, 500)
theo_fc = (1 + np.exp(-bohr_param))**-1

# Set up the figure canvas.
fig = plt.figure(figsize=(3.42, 5.55))
gs = gridspec.GridSpec(10, 1)
ax0 = fig.add_subplot(gs[:5])
ax1 = fig.add_subplot(gs[7:])
ax = [ax0, ax1]

# Format the various axes.
for a in ax:
    a.tick_params(labelsize=8)
ax[1].set_xlim([-0.5, len(dbls) - 0.5])
ax[1].set_xticks(np.arange(0, len(dbls)))
ax[1].set_xticklabels(dbls_labels)
ax[1].xaxis.grid(False)
ax[1].set_ylim([-5, 5])
ax[0].set_xlim([-10, 10])

# Add the appropriate labels.
ax[0].set_xlabel('Bohr parameter ($k_BT$)', fontsize=8)
ax[0].set_ylabel('fold-change', fontsize=8)
ax[1].set_ylabel(r'$\Delta F_{c \rightarrow \infty}$', fontsize=8)

# Plot the theoretical collapse curve.
_ = ax[0].plot(bohr_param, theo_fc, 'k-', lw=1.5, label='prediction')

# Add appropriate shading to the ΔF plot.
for i in range(len(dbls)):
    if i % 2 == 0:
        _ = ax[1].vlines(i, -5, 5, color='w', lw=35, alpha=0.75)

# Add the zero line.
_ = ax[1].hlines(0, -0.5, len(dbls) - 0.5, 'k',
                 linestyle=':', label='__nolegend__')

# Plot the data collapse.
grouped = data.groupby(['mutant', 'IPTGuM'])
for g, d in grouped:
    if 'Q294K' in g[0]:
        pass
    else:
        if g[0] == 'wt':
            face = 'w'
            dbl_mut = ['wt', 'wt']
        else:
            dbl_mut = g[0].split('-')
            face = colors[g[0]]
        edge = colors[g[0].upper()]

        if g[1] == 0.0:
            label = g[0].upper()
        else:
            label = '__nolegend__'

        # Compute the mean and error for the fold-change
        mean_fc = d['fold_change'].mean()
        err_fc = d['fold_change'].std() / np.sqrt(len(d))

        # Isoalte the parameters from the singles.
        ep_r = DNA_chains.iloc[DNA_ind]['{}_eps_r'.format(dbl_mut[0])]
        ka = IND_chains.iloc[IND_ind]['{}_Ka'.format(dbl_mut[1])]
        ki = IND_chains.iloc[IND_ind]['{}_Ki'.format(dbl_mut[1])]

        # Compute the bohr parameter.
        arch = mut.thermo.SimpleRepression(R=260, ep_r=ep_r, ka=ka, ki=ki,
                                           ep_ai=4.5, effector_conc=g[1]).bohr_parameter()

        # Plot the data collapse
        _ = ax[0].errorbar(arch, mean_fc, err_fc, color=colors[g[0].upper(
        )], markerfacecolor=face, markeredgecolor=edge, markeredgewidth=1, lw=1, markersize=6, label=label,
        fmt='.')
_ = ax[0].legend(loc='lower right', fontsize=8, handletextpad=0.05, handlelength=1)




# Compute the Delta F
def delta_F(mut_epr, mut_ka, mut_ki, wt_epr=-13.9, wt_ka=139, wt_ki=0.53, n=2,
            ep_ai=4.5):
    DNA_contrib = wt_epr - mut_epr
    numer = (1 + np.exp(-ep_ai) * (mut_ka / mut_ki)**n) 
    denom = (1 + np.exp(-ep_ai) * (wt_ka / wt_ki)**n)
    return DNA_contrib - np.log(numer / denom)

for i, m in enumerate(dbls):
    
    # Compute the measured value.
    meas_epr= DBL_chains.iloc[DBL_ind]['{}_eps'.format(m)]
    meas_ka= np.exp(-DBL_chains.iloc[DBL_ind]['{}_ka'.format(m)])
    meas_ki= np.exp(-DBL_chains.iloc[DBL_ind]['{}_ki'.format(m)])
    meas_df = delta_F(meas_epr, meas_ka, meas_ki)

    # Compute the predicted value.
    ep_r = DNA_chains.iloc[DNA_ind]['{}_eps_r'.format(m.split('-')[0])]
    ka = IND_chains.iloc[IND_ind]['{}_Ka'.format(m.split('-')[1])]
    ki = IND_chains.iloc[IND_ind]['{}_Ki'.format(m.split('-')[1])]
    pred_df = delta_F(ep_r, ka, ki)

    # Plot 'em
    ax[1].plot(i - 0.2, pred_df, 'o', ms=4, markeredgecolor=colors[m], markerfacecolor='w', label='__nolegend__')
    ax[1].plot(i + 0.2, meas_df, 'D', ms=4, markeredgecolor=colors[m], markerfacecolor='w', label='__nolegend__')

# Add the legends.
_ = ax[1].plot([], [], linestyle='none', marker='o', color='slategray', alpha=0.5, label='predicted')
_ = ax[1].plot([], [], linestyle='none', marker='D', color='slategray', alpha=0.5, label='measured')
_ = ax[1].legend(loc='upper left', fontsize=8)

plt.tight_layout()
plt.savefig('../../figures/fig5_vertical_plots.svg', bbox_inches='tight')
# %%
