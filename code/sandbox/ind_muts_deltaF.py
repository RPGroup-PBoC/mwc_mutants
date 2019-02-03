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

# Isolate the IND binding mutants and the wild-type.
IND = data[(data['class']=='IND') & (data['operator']=='O2')].copy()
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

IND_grouped = IND.groupby(['IPTGuM', 'mutant', 'repressors'])[['fold_change', 
                                    'delta_F']].agg(
                                    ('mean', 'sem')).reset_index()
# Set up the figure canvas. 
fig, ax = plt.subplots(2, 2, figsize=(6, 5))

# Plot the reference profile. 
for i in range(2):
    ax[i, 0].plot(c_range / Ka_0, ref_fc, '-', label='reference')
    ax[i, 1].plot(c_range / Ka_0, delF, '--', label='reference')

    ax[i, 0].errorbar(WT_grouped['IPTGuM'] / Ka_0,  WT_grouped['fold_change']['mean'], 
               WT_grouped['fold_change']['sem'], fmt='.', linestyle='none', capsize=1,
               lw=1, label='WT', ms=5)
    ax[i, 1].errorbar(WT_grouped['IPTGuM'] / Ka_0,  WT_grouped['delta_F']['mean'], 
               WT_grouped['delta_F']['sem'], fmt='.', linestyle='none', capsize=1,
               lw=1, label='WT', ms=5)

# Plot the summarized IND data. 
_colors = {m:c for m, c in zip(IND['mutant'].unique(), colors[2:])}
for g, d in IND_grouped.groupby(['mutant']):
    if g in ['F164T', 'Q294V']:
        row = 0
    else:
        row = 1 
    ax[row, 0].errorbar(d['IPTGuM'] / Ka_0,  d['fold_change']['mean'], 
               d['fold_change']['sem'], fmt='.', linestyle=':', capsize=1,
               lw=1, label=g,ms=5, color=_colors[g])
    ax[row,1].errorbar(d['IPTGuM'] / Ka_0,  d['delta_F']['mean'], 
        d['delta_F']['sem'], fmt='.', linestyle=':', capsize=1,
        lw=1, label='WT', ms=5, color=_colors[g])

# FOrmat and add labels       
for a in ax.ravel():
    a.set_xscale('log')
    a.set_xlabel('$c/ K_{A_0}$')

for i in range(2):
    ax[i, 0].set_ylabel('fold-change')
    ax[i, 1].set_ylabel('$\Delta F$ [$k_BT$]')
fig.text(0.0, 1, '(a)', fontsize=9)
fig.text(0.0, 0.5, '(b)', fontsize=9)

for i in range(2):
    ax[i, 0].set_ylim([-0.15, 1.15])
    ax[i, 1].set_ylim([-8, 8])
    ax[i, 0].legend(fontsize=8, handlelength=1)

plt.tight_layout()
plt.savefig('IND_mut_deltaF_r260.pdf', bbox_inches='tight')

