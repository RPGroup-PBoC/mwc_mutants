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
pboc_colors = mut.viz.color_selector('pboc')

# Load the data
data = pd.read_csv('../../data/csv/compiled_data.csv')
single_muts = data[(data['class'] != 'WT') & (data['class'] != 'DBL')]

# Restrict only to double mutants.
dbl = data[data['class'] == 'DBL'].copy()

# Generate a mean fold-change dataframe for easy plotting.
grouped = dbl.groupby(['mutant', 'IPTGuM']).apply(mut.stats.compute_mean_sem)
fc_df = pd.DataFrame(grouped).reset_index()

# Load the chains from the single mutant fits.
DNA_chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_DNA_strict.csv')
IND_chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_IND_strict.csv')

# Compute the stats for each.
DNA_idx = pd.DataFrame(DNA_chains.iloc[np.argmax(
    DNA_chains['lnprobability'])]).reset_index()
DNA_idx.columns = ['parameter', 'mode']
IND_idx = pd.DataFrame(IND_chains.iloc[np.argmax(
    IND_chains['lnprobability'])]).reset_index()
IND_idx.columns = ['parameter', 'mode']
modes = pd.concat([DNA_idx, IND_idx])

# Set some constants.
c_range = np.logspace(-8, -2, 500)
ka_wt = 139E-6
ki_wt = 0.53E-6
ep_R_wt = -13.9

#%%  Compute the Delta F.
dbl_muts = data[data['class']=='DBL']['mutant'].unique()
# Load the chains.
DNA_chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_DNA_strict.csv')
IND_chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_IND_strict.csv')
DBL_chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_DBL_strict.csv')
global_chains = pd.read_csv('../../data/mcmc/NB_emcee_mutants_global_strict.csv')

#%% Compute the statistics
DBL_stats = mut.stats.compute_statistics(DBL_chains, logprob_name='lnprobability')
global_stats = mut.stats.compute_statistics(global_chains, logprob_name='lnprobability')
stats = pd.concat([DNA_stats, IND_stats], ignore_index=True)
global_stats = pd.concat([global_stats, DBL_stats], ignore_index=True)

# %%
# Compute the WT bohr parameter.
bohr_c = np.logspace(-9, 9, 1000)
wt_eps_r = -13.9
wt_ka = 139E-6
wt_ki = 0.53E-6
wt_bohr = mut.thermo.SimpleRepression(R=260, ep_r=wt_eps_r, ka=wt_ka, ki=wt_ki, ep_ai=4.5,
                                       effector_conc=bohr_c).bohr_parameter()
pred_delta = {}
pred_c = {}
meas_delta = {}
meas_c = {}
# Loop through each mutant and calculated the predicted F.
for i, m in enumerate(dbl_muts):
    if ('Q294K' not in m) & ('Q294R' not in m):
        DNA_mut = m.split('-')[0]
        IND_mut = m.split('-')[1]
        if IND_mut != 'Q294K':
            pred_ep_r = stats[stats['parameter']=='{}_eps_r'.format(DNA_mut)]['mode'].values[0]
            pred_ka = np.exp(-stats[stats['parameter']=='{}_ka'.format(IND_mut)]['mode'].values[0])
            pred_ki = np.exp(-stats[stats['parameter']=='{}_ki'.format(IND_mut)]['mode'].values[0])
    # Compute the measured parameter
    meas_ep_r = global_stats[global_stats['parameter']=='{}_eps'.format(m)]['mode'].values[0]
    meas_ka = np.exp(-global_stats[global_stats['parameter']=='{}_ka'.format(m)]['mode'].values[0])
    meas_ki = np.exp(-global_stats[global_stats['parameter']=='{}_ki'.format(m)]['mode'].values[0])

    # Assemble the architectures and compute.  
    predicted_bohr = mut.thermo.SimpleRepression(R=260, ep_r=pred_ep_r, ka=pred_ka/1E6, ki=pred_ki/1E6, ep_ai=4.5,
                                              effector_conc=bohr_c).bohr_parameter()
    measured_bohr = mut.thermo.SimpleRepression(R=260, ep_r=meas_ep_r, ka=meas_ka/1E6, ki=meas_ki/1E6, ep_ai=4.5,
                                              effector_conc=bohr_c).bohr_parameter()
    pred_delta_c = wt_bohr - predicted_bohr
    pred_c[m] = pred_delta_c
    meas_delta_c = wt_bohr - measured_bohr
    meas_c[m] = meas_delta_c
    pred_ind = np.argmax(pred_delta_c)
    meas_ind = np.argmax(meas_delta_c)
    print('prediction', m,  bohr_c[pred_ind])
    print('measured', m, bohr_c[pred_ind])
    pred_delta[m] = pred_delta_c.sum()
    meas_delta[m] = meas_delta_c.sum()

# %% Set up the somewhat complicated figure.
colors = mut.viz.color_selector('mut')
fig = plt.figure(figsize=(3.5, 9))

# Set up the axes.
gs = gridspec.GridSpec(19, 3)
axes = [[], [], []]
for i in range(3):
    for j in range(0, 6, 2):
        _ax = fig.add_subplot(gs[j:j+2, i])
        _ax.tick_params(labelsize=6)
        _ax.set_xscale('log')
        _ax.set_xlim([1E-8, 1E-2])
        _ax.set_ylim([-0.1, 1.1])
        _ax.xaxis.set_ticks([1E-6, 1E-3])
        _ax.xaxis.set_tick_params()
        _ax.yaxis.set_ticks([0, 0.5, 1])
        axes[i].append(_ax)
ax_bohr = fig.add_subplot(gs[7:12, :])
ax_deltaF = fig.add_subplot(gs[16:, :])

# Properly label and format the axes.
dna_muts = ['Y20I', 'Q21A', 'Q21M']
for i in range(3):
    axes[i][0].set_title(
        dna_muts[i], backgroundcolor=pboc_colors['pale_yellow'], fontsize=8, y=1.04)
    axes[i][0].set_xticklabels([])
    axes[i][1].set_xticklabels([])
    if i > 0:
        axes[i][0].set_yticklabels([])
        axes[i][1].set_xticklabels([])
        axes[i][1].set_yticklabels([])
        axes[i][2].set_yticklabels([])

ind_muts = ['F164T', 'Q294K', 'Q294V']
for i in range(3):
    axes[0][i].text(-0.75, 0.65, ind_muts[i], fontsize=8, backgroundcolor=pboc_colors['pale_yellow'],
                    transform=axes[0][i].transAxes, rotation='vertical')
    if i == 1:
        axes[0][i].set_ylabel('fold-change', fontsize=8)
        axes[-1][i].set_xlabel('IPTG (M)', fontsize=8)
ax_bohr.tick_params(labelsize=8)
ax_bohr.set_xlabel('Bohr parameter ($k_BT$)', fontsize=8)
ax_bohr.set_ylabel('fold-change', fontsize=8)
ax_bohr.set_xlim([-10, 10])
ax_deltaF.tick_params(labelsize=8)
ax_deltaF.set_ylabel('$\Delta F$ $(k_BT)$', fontsize=8)

# Add panel labels.
fig.text(-0.05, 0.9, '(A)', fontsize=8)
fig.text(-0.05, 0.6, '(B)', fontsize=8)
fig.text(-0.05, 0.35, '(C)', fontsize=8)

# plot the double mutant data and predictions.
# Set up the axes locator.
mut_axes = {}
for i, dna in enumerate(dna_muts):
    for j, ind in enumerate(ind_muts):
        mut_axes['{}-{}'.format(dna, ind)] = axes[i][j]

#e Loop through the axes and plot the data theoretical predictions.
grouped = fc_df.groupby(['mutant'])
for g, d in grouped:
    axis = mut_axes[g]

    # oompute and plotthe theoretical prediction.
    if 'Q294K' not in g:
        # Get the stats.
        DNA_mut, IND_mut = g.split('-')
        ep_R = modes[modes['parameter'] ==
                     '{}_eps_r'.format(DNA_mut)]['mode'].values[0]
        ka = modes[modes['parameter'] ==
                   '{}_Ka'.format(IND_mut)]['mode'].values[0]
        ki = modes[modes['parameter'] ==
                   '{}_Ki'.format(IND_mut)]['mode'].values[0]

        # Set up the architecture.
        arch = mut.thermo.SimpleRepression(
            R=260, ep_r=ep_R, ep_ai=4.5, ka=ka/1E6, ki=ki/1E6, effector_conc=c_range)
        fc_theo = arch.fold_change()
        _ = axis.plot(c_range, fc_theo, lw=0.75, color=colors[g])

    _ = axis.errorbar(d['IPTGuM'] / 1E6, d['mean'], d['sem'],
                      ms=2, fmt='.', lw=1, color=colors[g])

# Plot the collapse curve.
# All data.
data = pd.read_csv('../../data/csv/compiled_data.csv')
for m in ['Q21A', 'Q21M', 'Y20I']:
    data = data[data['mutant'] != m + '-Q294K']
data = data[(data['class'] == 'DBL') | (data['class'] == 'WT')]
grouped = data.groupby(['mutant', 'repressors', 'IPTGuM'])
bohr = np.linspace(-10, 10, 500)
theo = (1 + np.exp(-bohr))**-1
ax_bohr.plot(bohr, theo, 'k-', label='prediction', zorder=1)

labeled = []
for g, d in grouped:
    split = g[0].split('-')
    if 'Q294K' not in split:
        if len(split) == 2:
            ep_R = modes[modes['parameter'] ==
                         '{}_eps_r'.format(split[0])]['mode'].values[0]
            ka = modes[modes['parameter'] == '{}_Ka'.format(
                split[1])]['mode'].values[0] * 1E-6
            ki = modes[modes['parameter'] == '{}_Ki'.format(
                split[1])]['mode'].values[0] * 1E-6
        else:
            if split[0] in dna_muts:
                ep_R = modes[modes['parameter'] ==
                             '{}_eps_r'.format(split[0])]['mode'].values[0]
                ka = ka_wt
                ki = ki_wt
            else:
                ep_R = ep_R_wt
                ka = modes[modes['parameter'] == '{}_Ka'.format(
                    split[0])]['mode'].values[0] * 1E-6
                ki = modes[modes['parameter'] == '{}_Ki'.format(
                    split[0])]['mode'].values[0] * 1E-6

        # Set up the architecture.
        arch = mut.thermo.SimpleRepression(R=g[1], ep_r=ep_R, ka=np.abs(
            ka), ki=np.abs(ki), ep_ai=4.5, effector_conc=g[2]/1E6)
        bohr_param = arch.bohr_parameter()

        # Compute the mean and SEM of the fold_change.
        mean_fc = d['fold_change'].mean()
        sem_fc = d['fold_change'].std() / np.sqrt(len(d))
        if split[0] == 'wt':
            face = 'w'
            edge = 'k'
        else:
            face = colors[g[0]]
            edge = colors[g[0]]
        if g[0] not in labeled:
            label = g[0].upper()
            labeled.append(g[0])
        else:
            label = '__nolegend__'
        _ = ax_bohr.errorbar(bohr_param, mean_fc, sem_fc, fmt='.', color=colors[g[0].upper(
        )], lw=1, linestyle='none', alpha=0.8, markerfacecolor=face, markeredgecolor=edge, markeredgewidth=1, markersize=7, label=label)
ax_bohr.legend(loc='upper left', handletextpad=0.1, fontsize=7.5, handlelength=0.75)

# Plot the delta F.
order = ['Y20I-F164T', 'Y20I-Q294V', 'Q21A-F164T', 'Q21A-Q294V', 'Q21M-F164T', 
        'Q21M-Q294V']
index = {o:i for i, o in enumerate(order)}
labels = ['Y20I\nF164T', 'Y20I\nQ294V', 'Q21A\nF164T', 'Q21A\nQ294V', 'Q21M\nF164T', 
        'Q21M\nQ294V']
# Define the colors 
ax = ax_deltaF
for i, delta in enumerate(pred_delta):
    mlabel = '__nolegend__'
    plabel = '__nolegend__'
    if ('Q294K' not in delta) & ('Q294R' not in delta):
        # Define the colors.
        color = colors[delta]
        _ = ax.plot(index[delta] + 0.2, pred_delta[delta], 'o', markeredgecolor=color, label=plabel, markerfacecolor='w', markersize=4, markeredgewidth=1.5)
        _ = ax.plot(index[delta] - 0.2, meas_delta[delta],  'D', markeredgecolor=color, label=mlabel, markerfacecolor='w', markersize=4, markeredgewidth=1.5)

# Format the axes further. 
ax.set_ylabel('$\int_{0}^\infty\Delta F$ dc', fontsize=8)
ax.set_xticks(np.arange(0, len(order), 1))
ax.set_xticklabels(labels)
ax.yaxis.grid(False)
for i in range(len(order)):
    if i%2 == 0:
        ax.vlines(i, -4500, 4500, color='w', lw=35, zorder=-1, alpha=0.75)
ax.set_ylim([-4500, 4500])
ax.set_xlim([-0.4, len(order) - 0.5])
_ = ax.plot([], [], 'o', color='slategray', label='predicted', alpha=0.5)
_ = ax.plot([], [], 'D', color='slategray', label='measured', alpha=0.5)
_ = ax.legend(loc='upper left', fontsize=8, handletextpad=0.03)
ax.hlines(0, -1, len(order)+1, color='k', linestyle=':')
plt.savefig('../../figures/doubles_data_collapse_vertical.svg', bbox_inches='tight')



#%%
fig, ax = plt.subplots(1, 1)
ax.set_xscale('log')
ax.set_ylabel('ΔF')
ax.set_xlabel('IPTG (M)')
for i, m in enumerate(pred_c):
    ax.plot(bohr_c, pred_c[m], color=colors[m], label=m)
    ax.plot(bohr_c, meas_c[m], color=colors[m], label='__nolegend__', linestyle=':')
ax.legend(loc='center right', fontsize=8)
plt.savefig('delta_f.pdf')
# plt.semilogx(bohr_c, wt_bohr, 'b')

#%%
DNA_stats
global_stats[global_stats['parameter']=='Q21M_eps']
global_stats[global_stats['parameter']=='Q21M-F164T_eps']