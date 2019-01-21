# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
constants = mut.thermo.load_constants()
colors = mut.viz.personal_style()

# Load the data. 
data = pd.read_csv('../../data/csv/compiled_data.csv')
DNA = data[((data['class']=='DNA') | (data['class']=='WT'))].copy()

# FOr empirical calculation of F, need to restrict to bounds of 0 and 1
# DNA.loc[(DNA['fold_change']> 1), 'fold_change'] = 1
# DNA.loc[(DNA['fold_change']< 0), 'fold_change'] = 0

c_range = np.logspace(-2, 4, 500)
c_range[0] = 0 

c_0 = 50

# Compute the measured F
DNA['measured_F'] = np.log((1 / DNA['fold_change']) - 1)

# Define the reference architecture
ref = mut.thermo.SimpleRepression(R= DNA['repressors'],
                                      ep_r=constants['O2'],
                                      ka=constants['Ka'],
                                      ki = constants['Ki'],
                                      ep_ai = constants['ep_AI'],
                                      effector_conc=c_0)
DNA['empirical_delta_F'] = -ref.bohr_parameter() - DNA['measured_F']

# Set up plotting things.
r, c = np.meshgrid(DNA['repressors'].unique(), c_range)
ref_arch = mut.thermo.SimpleRepression(R=r,
                                      ep_r=constants['O2'],
                                      ka=constants['Ka'],
                                      ki = constants['Ki'],
                                      ep_ai = constants['ep_AI'],
                                      effector_conc=c)

ref_arch_data = mut.thermo.SimpleRepression(R=DNA['repressors'],
                                      ep_r=constants['O2'],
                                      ka=constants['Ka'],
                                      ki = constants['Ki'],
                                      ep_ai = constants['ep_AI'],
                                      effector_conc=DNA['IPTGuM'])

DNA['predicted_delta_F'] = -ref.bohr_parameter() + ref_arch_data.bohr_parameter()
DNA['delta_delta_F'] = DNA['empirical_delta_F'] - DNA['predicted_delta_F']

fig, ax = plt.subplots(1, 3, figsize=(7, 2.5))

# Plot only the R=260 data. 

ax[0].plot(c[:, 1] / c_0 , ref_arch.fold_change()[:, 1])
ax[0].set_xlabel('$c / c_0$')
ax[0].set_ylabel('$c / c_0$')

grouped = DNA[(DNA['class']=='DNA') &
              (DNA['operator']=='O2')
             ].groupby(['mutant', 'repressors', 'IPTGuM'
                       ])['fold_change', 'delta_delta_F', 'predicted_delta_F',
                        ].agg(('mean', 'sem')).reset_index()


repressor_glyphs = {60:'D', 124:'X', 260:'o', 1220:'s'}
_colors = {'Q21M':colors[2], 'Y20I':colors[3], 'Q21A':colors[4], 'wt':colors[0]}

for g, d in grouped.groupby(['mutant', 'repressors']):
    if g[1] == 260:
        ax[0].errorbar(d['IPTGuM'] / c_0, d['fold_change']['mean'], d['fold_change']['sem'], fmt='o', 
                       color=_colors[g[0]], ms=4, lw=1, capsize=1, label=g[0])
        
# Plot the delta F
for g, d in grouped[grouped['repressors'] == 60].groupby('mutant'):
    ax[1].errorbar(d['IPTGuM'] / c_0, d['delta_delta_F']['mean'], 
               d['delta_delta_F']['sem'], marker='o', 
               color=_colors[g], ms=4, lw=1, linestyle='none', 
                capsize=1, label=g.upper())

ax[0].legend()
ax[0].set_ylim([-0.25, 1.25])
ax[1].set_ylim([-8, 8])
ax[0].set_xscale('symlog', linthreshx=c_range[2])
ax[1].set_xscale('symlog', linthreshx=c_range[2])


for g, d in DNA[(DNA['mutant']=='Y20I') & (DNA['repressors']==260)].groupby(['date', 'username']):
    plt.semilogx(d['IPTGuM'], d['fold_change'], label=g[0])
plt.legend()
                         