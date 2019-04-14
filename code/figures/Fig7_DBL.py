# -*- coding: utf-8 -*-
import sys
sys.path.insert(0,'../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.viz
import mut.thermo
import mut.stats
constants = mut.thermo.load_constants()
pboc = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load and prune data and deltaF
data = pd.read_csv('../../data/csv/summarized_data.csv')
data = data[data['class']=='DBL'].copy()
deltaF = pd.read_csv('../../data/csv/empirical_F_statistics.csv')
deltaF = deltaF[deltaF['class']=='DBL'].copy()


# Load the sampling information
epRA_samples = pd.read_csv('../../data/csv/DNA_binding_energy_samples.csv')
epRA_samps = epRA_samples[(epRA_samples['operator']=='O2') & 
                           (epRA_samples['repressors']==260)].copy()
kaki_epai_samples = pd.read_csv('../../data/csv/KaKi_epAI_samples.csv')
kaki_epai_samps = kaki_epai_samples[kaki_epai_samples['operator']=='O2'].copy()

# ##############################################################################
# FIGURE INSTANTIATION AND FORMATTING 
# ##############################################################################
fig, ax = plt.subplots(3, 7, figsize=(7.1, 2.8))
# Define the mutant axes
DNA_idx = {'Y20I':0, 'Q21A':1, 'Q21M':2}
IND_idx = {'F164T': 0, 'Q294V':1, 'Q294K':2}

for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)

# Disable middle column
for i in range(3):
    ax[i, 3].axis('off')

# Set scaling
for i in range(3):
    for j in range(3):
        ax[i, j].set_xscale('symlog', linthreshx=1E-3, linscalex=0.5)
        ax[i, j].set_xlim([1E-2, 1E4])
        ax[i, j].set_xlim([-0.001, 1E4])
        ax[i, j].set_xticks([0, 1E-1, 1E1, 1E3])
        ax[i, j].set_ylim([-0.1, 1.2])
        ax[i, j+4].set_ylim([-9, 9])
        ax[i, j].set_yticks([0, 0.5, 1])

# Set the titles and axis labels. 
for m, idx in IND_idx.items():
        ax[0, idx].set_title(m, fontsize=7, y=1.04, 
                           backgroundcolor=pboc['pale_yellow'])
        ax[0, idx + 4].set_title(m, fontsize=7, y=1.04, 
                           backgroundcolor=pboc['pale_yellow'])
        ax[-1, idx].set_xlabel('IPTG [µM]', fontsize=8)
        ax[-1, idx + 4].set_xlabel('$F^{(ref)}$ [$k_BT$]', fontsize=8)

# Remove unnecessary ticklabels
for i in range(2):
    ax[i,0].set_xticklabels([])
    ax[i,4].set_xticklabels([])
    ax[-1, i+1].set_yticklabels([])
    ax[-1, i+5].set_yticklabels([])
    for j in range(2):
        ax[j, i+1].set_xticklabels([])
        ax[j, i+1].set_yticklabels([])
        ax[j, i+5].set_xticklabels([])
        ax[j, i+5].set_yticklabels([])

# Mutant Identifiers
for m, idx in DNA_idx.items():
        ax[idx, 0].text(-0.85, 0.62, m, fontsize=7,rotation='vertical', 
            backgroundcolor=pboc['pale_yellow'], transform=ax[idx, 0].transAxes)
        ax[idx, 4].text(-0.85, 0.62, m, fontsize=7,rotation='vertical', 
            backgroundcolor=pboc['pale_yellow'], transform=ax[idx, 4].transAxes)
        ax[idx, 0].set_ylabel('fold-change', fontsize=7, labelpad=0.1)
        ax[idx, 4].set_ylabel('$\Delta F$ [$k_BT$]', fontsize=7, labelpad=0.01)

# Panel Labels
fig.text(0, 0.94, '(A)', fontsize=8)
fig.text(0.48, 0.94, '(B)', fontsize=8)
# ##############################################################################
# INDUCTION DATA 
# ##############################################################################
for dna, dna_idx in DNA_idx.items():
    for ind, ind_idx in IND_idx.items():
        # Get the mutant from the summarized data. 
        _data = data[data['mutant']==f'{dna}-{ind}']

        ax[dna_idx, ind_idx].errorbar(_data['IPTGuM'], _data['mean'], _data['sem'],
            fmt='.', lw=1, capsize=1, linestyle='none', color=pboc['purple'],
            ms=3)

# ##############################################################################
# DELTA F DATA
# ##############################################################################
# compute the reference
arch = mut.thermo.SimpleRepression(R=260, ep_r=-13.9, ka=153, ki=0.53, n_sites=2,
    effector_conc=data['IPTGuM'].unique(), ep_ai=4.5)
ref = arch.bohr_parameter()
for dna, dna_idx in DNA_idx.items():
    for ind, ind_idx in IND_idx.items():
        _data = deltaF[(deltaF['mutant']==f'{dna}-{ind}') & 
                       (deltaF['parameter']=='delta_bohr_corrected')]

        ax[dna_idx, ind_idx + 4].plot(ref, _data['median'], '.', 
                                color=pboc['purple'], ms=3)
        ax[dna_idx, ind_idx + 4].vlines(ref, _data['hpd_min'], _data['hpd_max'],
                                        lw=1, color=pboc['purple'])

# ##############################################################################
# THEORY CURVES
# ##############################################################################
c_range = np.logspace(-3, 4, 200)
c_range[0] = 0
ref_bohr = mut.thermo.SimpleRepression(R=260, ep_r=-13.9, ka=139, ki=0.53,
                                        n_sites=2, effector_conc=c_range,
                                        ep_ai=4.5).bohr_parameter()
n_draws = int(1E4)
for dna, dna_idx in DNA_idx.items():
    for ind, ind_idx in IND_idx.items():
        _allo_samps = kaki_epai_samps[kaki_epai_samps['mutant']==ind]
        _epRA_samps = epRA_samps[epRA_samps['mutant']==dna]
        epRA_draws = np.random.choice(_epRA_samps['ep_RA'].values, replace=True,
                                      size=n_draws)
        allo_draws = _allo_samps.sample(n=n_draws, replace=True)
        fc_cred_region = np.zeros((2, len(c_range)))
        bohr_cred_region = np.zeros((2, len(c_range)))
        for i, c in enumerate(c_range):
            arch = mut.thermo.SimpleRepression(R=260, ep_r=epRA_draws,
                            ka=allo_draws['Ka'], ki=allo_draws['Ki'], 
                            ep_ai=allo_draws['ep_AI'], n_sites=2, effector_conc=c)
            fc_cred_region[:, i] = mut.stats.compute_hpd(arch.fold_change(), 0.95)
            bohr_cred_region[:, i] = mut.stats.compute_hpd(ref_bohr[i] -\
                                                     arch.bohr_parameter(), 0.95)

        ax[dna_idx, ind_idx].fill_between(c_range, fc_cred_region[0, :], 
                                        fc_cred_region[1, :], alpha=0.5, 
                                        color=pboc['purple'])
        ax[dna_idx, ind_idx + 4].fill_between(ref_bohr, bohr_cred_region[0, :], 
                                        bohr_cred_region[1, :], alpha=0.5, 
                                        color=pboc['purple'])

plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.savefig('../../figures/Fig7_DBL_deltaF.pdf', bbox_inches='tight')