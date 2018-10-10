# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.thermo
import mut.viz
import mut.stats
pboc = mut.viz.color_selector('pboc')
colors = mut.viz.color_selector('mut')
mut.viz.plotting_style()
constants = mut.thermo.load_constants()


# Load the summarized and sampling data
data = pd.read_csv('../../data/csv/summarized_data.csv')
data = data[data['class']=='DBL'].copy()
DBL_fit = pd.read_csv('../../data/csv/Fig5_O2_DBL_samples.csv')
DBL_stats = pd.read_csv('../../data/csv/Fig5_O2_DBL_stats.csv')
DBL_pred = pd.read_csv('../../data/csv/Fig5_O2_DBL_predicted_Bohr.csv')
# DBL_fit['mutant'] = DBL_fit['DNA_mutant'] + '-' + DBL_fit['IND_mutant']

# Compute the standard collapse curve.
bohr_range = np.linspace(-8, 12, 200)
collapse = (1 + np.exp(-bohr_range))**-1

# Compute the wt bohr.
wt_bohr = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=constants['O2'],
                                     ka=constants['Ka'], ki=constants['Ki'],
                                     ep_ai=constants['ep_AI'], n_sites=constants['n_sites'],
                                     n_ns=constants['Nns'], effector_conc=data['IPTGuM'].unique()).bohr_parameter()

# Instantiate the figure, 
fig = plt.figure(figsize=(3.42, 6.5))
gs = gridspec.GridSpec(6, 1)
ax0 = fig.add_subplot(gs[:3, 0])
ax1 = fig.add_subplot(gs[5:, 0])

# Create the positional indices
pos = {m:i for i, m in enumerate(data['mutant'].unique())}

# Add appropriate labels and formatting.
for a in [ax0, ax1]:
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)

ax1.set_xticks(list(pos.values()))
ax1.set_xticklabels(pos.keys(), rotation=90)
ax1.set_ylim([-1, 40])
ax0.set_xlim([-8, 12])
ax0.set_xlabel('Bohr parameter [$k_BT$]', fontsize=8)
ax0.set_ylabel('fold-change', fontsize=8)
ax1.set_ylabel(r'$\chi^2_{\Delta F}$', fontsize=8)

# Plot the collapse curve. 
ax0.plot(bohr_range, collapse, 'k-', label='__nolegend__')

# Plot the measured bohr parameter
for g, d in DBL_pred.groupby('mutant'):
    d = d.sort_values('IPTGuM')
    fc = data[data['mutant']==g].sort_values('IPTGuM')
    ax0.errorbar(d['mode'], fc['mean'], fc['sem'], linestyle='none',
                lw=1, fmt='o', color=colors[g], label=g, ms=3)
    ax0.hlines(fc['mean'], d['hpd_min'], d['hpd_max'], lw=1, color=colors[g], label='__nolegend__')


for g, d in DBL_pred.groupby('mutant'):
    d = d.sort_values(['IPTGuM'])
    ssq_mode = np.sum((d['mode'].values - wt_bohr)**2) / len(d)
    ssq_max = np.sum((d['hpd_max'].values - wt_bohr)**2) / len(d)
    ssq_min = np.sum((d['hpd_min'].values - wt_bohr)**2) / len(d)

    _ = ax1.plot(pos[g] - 0.17, ssq_mode, 'o', markerfacecolor='w',
                markeredgecolor=colors[g], markeredgewidth=1, markersize=5,
                label='__nolegend__')
    _ = ax1.vlines(pos[g] - 0.17, ssq_min, ssq_max, lw=1, color=colors[g],
                  label='__nolegend__')

# Plot the measured sum square residual.
for i, m in enumerate(data['mutant'].unique()):
    # Compute the bohr parameter.
    _bohr_mode = mut.thermo.SimpleRepression(R=constants['RBS1027'],
                                    ep_r=DBL_stats[DBL_stats['parameter']=='ep_RA.{}'.format(m)]['mode'].values[0],
                                    ka=DBL_stats[DBL_stats['parameter']=='Ka.{}'.format(m)]['mode'].values[0],
                                    ki=DBL_stats[DBL_stats['parameter']=='Ki.{}'.format(m)]['mode'].values[0],
                                    ep_ai=constants['ep_AI'], effector_conc=data['IPTGuM'].unique(),
                                    n_sites=constants['n_sites'], n_ns=constants['Nns']).bohr_parameter()
    _bohr_min = []
    _bohr_max = []
    for c in np.sort(data['IPTGuM'].unique()):
        bohr = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=DBL_fit['ep_RA.{}'.format(m)],
                                          ka=DBL_fit['Ka.{}'.format(m)], ki=DBL_fit['Ki.{}'.format(m)],
                                          ep_ai=constants['ep_AI'], n_sites=constants['n_sites'],
                                          n_ns=constants['Nns'], effector_conc=c).bohr_parameter()
        low, high = mut.stats.compute_hpd(bohr, 0.95)
        _bohr_min.append(low)
        _bohr_max.append(high)

    ssq_median = np.sum((_bohr_mode - wt_bohr)**2) / len(_bohr_mode)
    ssq_max = np.sum((_bohr_max - wt_bohr)**2) / len(_bohr_mode)
    ssq_min = np.sum((_bohr_min - wt_bohr)**2) / len(_bohr_mode)
    _ = ax1.plot(pos[m] + 0.17, ssq_median, 'D', markerfacecolor='w',
                markeredgecolor=colors[m], markeredgewidth=1, markersize=5, label='__nolegend__')
    _ = ax1.vlines(pos[m] + 0.17, ssq_min, ssq_max, lw=1, color=colors[m], label='__nolegend__')

for m, i in pos.items():
    if i%2 == 0:
        _ = ax1.vlines(i, -1, 50, color=pboc['pale_yellow'], lw=25, zorder=-1)
ax0.legend(fontsize=7, handletextpad=0.1)
ax1.plot([], [], 'o', markerfacecolor='w', markeredgecolor='k', markersize=5, markeredgewidth=1,
        label='predicted')
ax1.plot([], [], 'D', markerfacecolor='w', markeredgecolor='k', markersize=5, markeredgewidth=1,
        label='measured')
ax1.legend(loc='upper left', fontsize=7, ncol=2)
plt.savefig('Fig5_collapse.svg', bbox_inches='tight')
