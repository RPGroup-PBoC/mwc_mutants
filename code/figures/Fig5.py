# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.thermo
import mut.stats
import mut.viz
color = mut.viz.color_selector('mut')
pboc = mut.viz.color_selector('pboc')
mut.viz.plotting_style()
constants = mut.thermo.load_constants()

# Load the data. 
indiv_samps = pd.read_csv('../../data/csv/Fig2_O2_DBL_samples.csv')
indiv_stats = pd.read_csv('../../data/csv/Fig2_O2_DBL_stats.csv')
predicted = pd.read_csv('../../data/csv/Fig5_O2_DBL_predicted_deltaBohr.csv')o

# Load the various statistics
epRA_stats = pd.read_csv('../../data/csv/Fig2_O2_DNA_binding_energy_stats.csv')
kaki_stats = pd.read_csv('../../data/csv/Fig3_O2_KaKi_only_stats.csv')
allo_stats = pd.read_csv('../../data/csv/Fig3_O2_KaKi_epAI_stats.csv')

# Instantiate the figure
fig, ax = plt.subplots(1, 2, figsize=(7,4))

# Assign the labels and format. 
for a in ax:
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)

ax[0].set_ylabel('fold-change', fontsize=8)
ax[0].set_xlabel('Bohr parameter [$k_BT$]', fontsize=8)
ax[1].set_ylabel(r'$\Delta F_{c\rightarrow\infty}$', fontsize=8)


# Compute and plot the master curve
bohr_range = np.linspace(-8, 8, 200)
theo = (1 + np.exp(-bohr_range))**-1
ax[0].plot(bohr_range, theo, 'k-', label='__nolegend__')

# For each mutant, compute the most-likely value of bohr parameter 
 
for g, d in 
predicted
