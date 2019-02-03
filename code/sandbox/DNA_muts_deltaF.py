# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
constants = mut.thermo.load_constants()
colors = mut.viz.personal_style()

# Load the data. 
data = pd.read_csv('../../data/csv/compiled_data.csv')
data.dropna(inplace=True)

# Compute the reference strain bohr. 
epRA_0 = np.array([constants[op] for op in data['operator'].values])
Ka_0 = constants['Ka']
Ki_0 = constants['Ki']
epAI_0 = constants['ep_AI']
R_0 = data['repressors']
n_sites = 2
c_0 = data['IPTGuM']
ref_arch = mut.thermo.SimpleRepression(R=R_0, effector_conc=c_0, ka=Ka_0, ki=Ki_0, 
                                   ep_ai=epAI_0, ep_r=epRA_0, n_sites=n_sites) 
data['reference_F'] = -ref_arch.bohr_parameter()

# Calculate the empirical F and delta F
data['empirical_F'] = np.log((1 / data['fold_change']) - 1)
data['delta_F'] = data['reference_F'] - data['empirical_F']

# Isolate the DNA binding mutants and the wild-type.
DNA = data[(data['class']=='DNA') & (data['operator']=='O2')].copy()
WT = data[(data['class']=='WT') & (data['operator']=='O2')].copy()

# Compute the reference strain for plotting purposes
c_range = np.logspace(-2, 4, 500)
ref_arch = mut.thermo.SimpleRepression(R=260, effector_conc=c_range, ka=Ka_0, ki=Ki_0, 
                                   ep_ai=epAI_0, ep_r=constants['O2'], n_sites=n_sites) 
ref_fc = ref_arch.fold_change()
delF = np.zeros_like(c_range)

# Compute the summarized data. 
WT_grouped = WT.groupby(['IPTGuM'])[['fold_change', 
                                    'delta_F']].agg(
                                    ('mean', 'sem')).reset_index()

DNA_grouped = DNA.groupby(['IPTGuM', 'mutant', 'repressors'])[['fold_change', 
                                    'delta_F']].agg(
                                    ('mean', 'sem')).reset_index()

# Set up the figure canvas. 
fig, ax = plt.subplots(1, 2, figsize=(6, 3))

# Plot the reference profile. 
ax[0].plot(c_range / Ka_0, ref_fc, '-', label='reference')
ax[1].plot(c_range / Ka_0, delF, '--', label='reference')

# Plot the summarized wt data. 
ax[0].errorbar(WT_grouped['IPTGuM'] / Ka_0,  WT_grouped['fold_change']['mean'], 
               WT_grouped['fold_change']['sem'], fmt='.', linestyle='none', capsize=1,
               lw=1, label='WT')
ax[1].errorbar(WT_grouped['IPTGuM'] / Ka_0,  WT_grouped['delta_F']['mean'], 
               WT_grouped['delta_F']['sem'], fmt='.', linestyle='none', capsize=1,
               lw=1, label='WT')

# Plot the summarized DNA data. 
glyphs = ['s', 'D', 'o', 'X']
reps = {r:g for r, g in zip(np.sort(data['repressors'].unique()), glyphs)}
for g, d in DNA_grouped.groupby(['mutant', 'repressors']):
    if g[1] == 260:
        ax[0].errorbar(d['IPTGuM'] / Ka_0,  d['fold_change']['mean'], 
               d['fold_change']['sem'], fmt='.', linestyle=':', capsize=1,
               lw=1, label=g[0])
        ax[1].errorbar(d['IPTGuM'] / Ka_0,  d['delta_F']['mean'], 
        d['delta_F']['sem'], fmt=reps[g[1]], linestyle=':', capsize=1,
        lw=1, label='WT', ms=3)

# FOrmat and add labels       
for a in ax:
    a.set_xscale('log')
    a.set_xlabel('$c/ K_{A_0}$')
ax[0].set_ylabel('fold-change')
ax[1].set_ylabel('$\Delta F$ [$k_BT$]')
fig.text(0.0, 1, '(a)', fontsize=9)
fig.text(0.5, 1, '(b)', fontsize=9)

ax[0].set_ylim([-0.15, 1.15])
ax[1].set_ylim([-8, 8])
ax[0].legend(fontsize=8, handlelength=1)

plt.tight_layout()
plt.savefig('dna_mut_deltaF_r260.pdf', bbox_inches='tight')

