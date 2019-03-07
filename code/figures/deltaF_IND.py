# -*- coding: utf-8 -*-
import sys 
sys.path.insert(0, '../../')
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import mut.thermo
import mut.viz 
mut.viz.plotting_style()
colors = mut.viz.color_selector('mut')
pboc = mut.viz.color_selector('pboc')
constants = mut.thermo.load_constants()

# Load the data and restrict to the inducer binding mutants. 
data = pd.read_csv('../../data/csv/empirical_F_statistics.csv')
IND = data[data['class']=='IND'].copy()

# Load the statistics and sampling information. 
kaki_only_stats = pd.read_csv('../../data/csv/KaKi_only_summary.csv')
kaki_only_stats = kaki_only_stats[kaki_only_stats['operator']=='O2'].copy()
kaki_only_samples= pd.read_csv('../../data/csv/KaKi_only_samples.csv')
kaki_only_samples = kaki_only_samples[kaki_only_samples['operator']=='O2'].copy()
kaki_epAI_stats = pd.read_csv('../../data/csv/KaKi_epAI_summary.csv')
kaki_epAI_stats = kaki_epAI_stats[kaki_epAI_stats['operator']=='O2'].copy()
kaki_epAI_samples= pd.read_csv('../../data/csv/KaKi_epAI_samples.csv')
kaki_epAI_samples = kaki_epAI_samples[kaki_epAI_samples['operator']=='O2'].copy()

# Compute the wild-type bohr parameter and pact
ep_r = [constants[op] for op in IND['operator']]
wt_bohr = mut.thermo.SimpleRepression(R=260, ep_r=ep_r,
                                      ka=constants['Ka'], ki=constants['Ki'],
                                      ep_ai=constants['ep_AI'], 
                                  effector_conc=IND['IPTGuM']).bohr_parameter()
c_range = np.logspace(-2, 4, 200)
c_range[0] = 0
wt_pact = mut.thermo.MWC(ka=constants['Ka'], ki=constants['Ki'],
                         ep_ai=constants['ep_AI'], 
                         effector_conc=c_range).pact()

# Instantiate the figure. 
fig, ax = plt.subplots(3, 4, figsize=(7, 4.5), sharex=True, sharey=True)
mut_idx = {'F164T': 0, 'Q294V': 1, 'Q294K': 2, 'Q294R': 3}
op_idx = {'O1': 0, 'O2': 1, 'O3': 2}

# ####################################
# DATA
# ####################################
op_glyphs = {'O1': '^', 'O2':'o', 'O3': 's'}
for g, d in IND.groupby(['mutant', 'operator']):
    _d = d[d['parameter']=='delta_bohr_corrected']
    _ax = ax[op_idx[g[1]], mut_idx[g[0]]]
    if g[1] == 'O2':
        face = 'w'
    else:
        face = colors[g[0]]
    _ax.plot(_d['IPTGuM'], _d['mean'], color=colors[g[0]], markerfacecolor=face,
            marker=op_glyphs[g[1]], ms=4, linestyle='none')
    _ax.vlines(_d['IPTGuM'], _d['hpd_min'], _d['hpd_max'], lw=0.75, color=colors[g[0]])

# ###################################
# THEORY
# ######################################
# epistasis free case
for d in IND['mutant'].unique():
    # Isolate the relevant samples and statistics
    if d in ['F164T', 'Q294V']:
        _samples = kaki_only_samples[kaki_only_samples['mutant']==d]
        _samples['ep_AI'] = constants['ep_AI']
        median = kaki_only_stats[kaki_only_stats['mutant']==d][['parameter', 'median']]
        Ka = median[median['parameter']=='Ka']['median'].values[0]
        Ki = median[median['parameter']=='Ki']['median'].values[0]
        ep_AI = constants['ep_AI']
    else:
        _samples = kaki_epAI_samples[kaki_epAI_samples['mutant']==d]
        median = kaki_epAI_stats[kaki_epAI_stats['mutant']==d][['parameter', 'median']]
        Ka = median[median['parameter']=='Ka']['median'].values[0]
        Ki = median[median['parameter']=='Ki']['median'].values[0]
        ep_AI = median[median['parameter']=='ep_AI']['median'].values[0]


    for i in range(3):
        _ax = ax[i, mut_idx[d]]
        # Plot the credible regions
        cred_region = np.zeros((2, len(c_range)))
        for j, c in enumerate(c_range):
            pact = mut.thermo.MWC(ka=_samples['Ka'], ki=_samples['Ki'],
                                      ep_ai=_samples['ep_AI'],
                                      effector_conc=c).pact()
            deltaF = -np.log(wt_pact[j] / pact)
            cred_region[:, j] = mut.stats.compute_hpd(deltaF, 0.95)
        _ = _ax.fill_between(c_range, cred_region[0, :], cred_region[1,:],
            color=colors[d], alpha=0.4, zorder=1)

        # Plot the medians
        pact = mut.thermo.MWC(ka=Ka, ki=Ki, ep_ai=ep_AI,
                              effector_conc=c_range).pact()
        deltaF = -np.log(wt_pact/pact)
        _ = _ax.plot(c_range, deltaF, lw=0.75, color=colors[d], zorder=2)


for a in ax.ravel():
    a.hlines(0, -1, 1E4, lw=0.75, linestyle='--', color='k', zorder=1)

# #####################################
# LABELING
# #####################################
for i in range(3):
    ax[i, 0].set_ylabel('$\Delta F$ [$k_BT$]', fontsize=8)
for i in range(4):
    ax[2, i].set_xlabel('IPTG [µM]', fontsize=8)
for d, ind in mut_idx.items():
    ax[0, ind].set_title(d, fontsize=9, backgroundcolor=pboc['pale_yellow'],
                        y=1.04)
for op, ind in op_idx.items():
    ax[ind, 0].text(-0.58, 0.47, f'{op}',
        rotation='vertical', transform=ax[ind, 0].transAxes, fontsize=9,
        backgroundcolor=pboc['pale_yellow'])
fig.text(-0.03, 0.6, 'operator identity', fontsize=10, rotation='vertical',
        backgroundcolor='#E3DBD0')
fig.text(0.4, 0.99, 'inducer binding mutant', fontsize=10, 
        backgroundcolor='#E3DBD0')
# ####################################
# FORMATTING
# ####################################
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=9)
    a.yaxis.set_tick_params(labelsize=9)
    a.set_xscale('symlog')
    a.set_xlim([-1, 1E4])
    a.set_ylim([-8, 8])
plt.subplots_adjust(wspace=0.08, hspace=0.08)
plt.savefig('../../figures/Fig5_IND_deltaF.pdf', bbox_inches='tight')
