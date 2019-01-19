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
IND = data[(data['class']=='IND') | (data['class']=='WT')].copy()
c_range = np.logspace(-2, 4, 500)
c_0 = 50

# Compute the measured F
IND['measured_F'] = np.log((1 / IND['fold_change']) - 1)

# Compute the empirical delta F
# Define the reference architecture
ref = mut.thermo.SimpleRepression(R=constants['RBS1027'],
                                      ep_r=constants['O2'],
                                      ka=constants['Ka'],
                                      ki = constants['Ki'],
                                      ep_ai = constants['ep_AI'],
                                      effector_conc=c_0)
IND['empirical_delta_F'] = -ref.bohr_parameter() - IND['measured_F']

# Set up plotting things.
ref_arch = mut.thermo.SimpleRepression(R=constants['RBS1027'],
                                      ep_r=constants['O2'],
                                      ka=constants['Ka'],
                                      ki = constants['Ki'],
                                      ep_ai = constants['ep_AI'],
                                      effector_conc=c_range)

ref_arch_data = mut.thermo.SimpleRepression(R=constants['RBS1027'],
                                      ep_r=constants['O2'],
                                      ka=constants['Ka'],
                                      ki = constants['Ki'],
                                      ep_ai = constants['ep_AI'],
                                      effector_conc=IND['IPTGuM'])
IND['predicted_delta_F'] = -ref.bohr_parameter() + ref_arch_data.bohr_parameter()
fig, ax = plt.subplots(4, 3, figsize=(7, 7))
for i in range(4):
    ax[i, 0].set_xscale('log')
    ax[i, 1].set_xscale('log')
    ax[i, 0].plot(c_range / c_0, ref_arch.fold_change(), 'k-', lw=1)
    ax[i, 1].plot(c_range / c_0, -ref.bohr_parameter() + ref_arch.bohr_parameter(), '-', color='rebeccapurple', lw=1)
    ax[i, 2].plot([-8, 8], [-8, 8], 'k:', lw=1)


# Plot the empirical deltaF
axes = {'F164T':0, 'Q294V':1, 'Q294K':2, 'Q294R':3}
grouped  = IND[(IND['class']=='IND') &
               (IND['repressors']==260) &
               (IND['operator']=='O2')
              ].groupby(['IPTGuM', 'mutant'])['empirical_delta_F', 
                                              'predicted_delta_F',
                                             'fold_change'].agg(('mean', 'sem')).reset_index()

    
for g, d in grouped.groupby('mutant'):
    ax[axes[g], 0].errorbar(d['IPTGuM'] / c_0, d['fold_change']['mean'], d['fold_change']['sem'],
                           lw=1, color=colors[3], capsize=1, marker='.', markersize=3,
                           linestyle='none')
    ax[axes[g], 1].errorbar(d['IPTGuM'] / c_0, d['empirical_delta_F']['mean'], d['empirical_delta_F']['sem'],
                           lw=1, color=colors[3], capsize=1, marker='.', markersize=3,
                           linestyle='none')

    ax[axes[g], 2].errorbar(d['predicted_delta_F']['mean'], d['empirical_delta_F']['mean'], 
                            d['empirical_delta_F']['sem'],
                           lw=1, color=colors[3], capsize=1, marker='.', markersize=3,
                           linestyle='none')       