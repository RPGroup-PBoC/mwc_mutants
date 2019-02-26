# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
from tqdm import tqdm
np.random.seed(666)
mut.viz.plotting_style()
colors = mut.viz.color_selector('mut')
pboc = mut.viz.color_selector('pboc')
constants = mut.thermo.load_constants()

# Load the data and restrict to the double mutants.
data = pd.read_csv('../../data/csv/empirical_F_statistics.csv')
DBL = data[data['class']=='DBL'].copy()
DBL.drop_duplicates(inplace=True)

# Load the necssary statistics. 
epRA_samples = pd.read_csv('../../data/csv/DNA_binding_energy_samples.csv')
epRA_samples = epRA_samples[(epRA_samples['operator']=='O2') & 
                            (epRA_samples['repressors']==260)]
kaki_only_samples = pd.read_csv('../../data/csv/KaKi_only_samples.csv')
kaki_only_samples = kaki_only_samples[kaki_only_samples['operator']=='O2']
kaki_only_samples['ep_AI'] = constants['ep_AI']
kaki_epAI_samples = pd.read_csv('../../data/csv/KaKi_epAI_samples.csv')
kaki_epAI_samples = kaki_epAI_samples[kaki_epAI_samples['operator']=='O2']

# Draw from distributions and compute credible regions for each double mut. 
n_draws = int(1E4)
c_range = DBL['IPTGuM'].unique()
c_range[0] = 0
draw_dfs = []
for d in tqdm(DBL['mutant'].unique()):
    dna_mut, ind_mut = d.split('-')
    # Draw the samples. 
    epRA_draws = epRA_samples[
        epRA_samples['mutant']==dna_mut]['ep_RA'].sample(n_draws, replace=True)
    if ind_mut == 'Q294K':
        _kaki_samples = kaki_epAI_samples
    else:
        _kaki_samples = kaki_only_samples
    ind_draws = _kaki_samples[_kaki_samples['mutant']==ind_mut][
                ['Ka', 'Ki', 'ep_AI']].sample(n_draws, replace=True)

    # Compute the credible regions. 
    cred_region = np.zeros((2, len(c_range)))
    median =np.zeros(len(c_range))
    for i, c in enumerate(c_range):
        mut_arch = -mut.thermo.SimpleRepression(R=260, ep_r=epRA_draws.values,
                                    ka=ind_draws['Ka'], ki=ind_draws['Ki'],
                                    ep_ai=ind_draws['ep_AI'],
                                    effector_conc=c).bohr_parameter()
        cred_region[:, i] = mut.stats.compute_hpd(mut_arch, 0.95)
        median[i] = np.median(mut_arch)

    # Assemble the data frame.
    _df = pd.DataFrame(np.array([c_range, median, cred_region[0, :], cred_region[1, :]]).T,
                       columns=['IPTGuM', 'median', 'bohr_min', 'bohr_max'])
    _df['mutant'] = d
    draw_dfs.append(_df)
draw_df = pd.concat(draw_dfs)
ind_draws['ep_AI']
j
# Instantiate the figure. 
fig, ax = plt.subplots(3, 3, figsize=(6, 6), sharex=True, sharey=True)
DNA_idx = {'Q21M':0, 'Q21A':1, 'Y20I':2}
IND_idx = {'F164T':0, 'Q294V': 1, 'Q294K': 2}

# ####################################
# DATA
# #####################################

wt_bohr = -mut.thermo.SimpleRepression(R=260, ep_r=constants['O2'],
                                      ka=constants['Ka'], ki=constants['Ki'],
                                      ep_ai=constants['ep_AI'], 
                                  effector_conc=c_range).bohr_parameter()

for g, d in DBL.groupby(['mutant']):
    mutant = draw_df[draw_df['mutant']==g]
    mutant.sort_values(by='IPTGuM', inplace=True)
    d.sort_values(by='IPTGuM', inplace=True)
    d['delta_F'] = mutant['median'].values - d['bohr_median']
    d['delta_F_min'] = mutant['median'].values - d['bohr_min']
    d['delta_F_max'] = mutant['median'].values - d['bohr_max']
    dna_mut, ind_mut = g.split('-')
    _ax = ax[DNA_idx[dna_mut], IND_idx[ind_mut]]

    # Plot the data and error bars
    _ = _ax.plot(d['IPTGuM'], d['delta_F'], 'o', color=colors[g], ms=4)
    _ = _ax.vlines(d['IPTGuM'], d['delta_F_min'], d['delta_F_max'], lw=0.75,
                color=colors[g])

# #####################################
# LABELING
# #####################################
for i in range(3):
    ax[i, 0].set_ylabel('$\Delta F$ [$k_BT$]', fontsize=9)
    ax[-1, i].set_xlabel('IPTG [µM]', fontsize=9)
for ind, i in IND_idx.items():
    ax[0, i].set_title(ind, fontsize=9, backgroundcolor=pboc['pale_yellow'],
                    y=1.02)
for dna, i in DNA_idx.items():
    ax[i, 0].text(-0.57, 0.54, dna, fontsize=9, 
                  backgroundcolor=pboc['pale_yellow'], rotation='vertical',
                  transform=ax[i, 0].transAxes)


# #####################################
# FORMATTING
# #####################################
for a in ax.ravel():
    a.set_xscale('symlog')
    a.hlines(0, -1, 1E4, color='k', lw=0.75, linestyle='--')
    a.set_ylim([-9, 10])
    a.set_xlim([-0.2, 1E4])
    a.xaxis.set_tick_params(labelsize=9)
    a.yaxis.set_tick_params(labelsize=9)
plt.subplots_adjust(wspace=0.08, hspace=0.12)
fig.text(-0.07, 0.6, 'DNA binding mutation', fontsize=10, rotation='vertical',
backgroundcolor='#E3DBD0')
fig.text(0.35, 0.95, 'inducer binding mutation', fontsize=10,
backgroundcolor='#E3DBD0')
plt.savefig('../../data/csv/Fig8_DBL_deltaF_predref.pdf', bbox_inches='tight')