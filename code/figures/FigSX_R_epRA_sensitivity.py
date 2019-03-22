# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
import imp
imp.reload(mut.viz)
constants = mut.thermo.load_constants()
pboc = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Define the ranges
R_range = np.logspace(1, 3, 500)
epRA_range = np.linspace(-20, -5, 500)

# Set up the figure
fig, ax = plt.subplots(1, 2, figsize=(4.5, 2))
_ = ax[0].semilogx(R_range, 1/R_range, lw=2, color=pboc['red'])
_ = ax[1].plot(epRA_range, np.ones_like(epRA_range) * -1, lw=2, 
            color=pboc['red'])

# Add labels 
_ = ax[0].set_xlabel('repressors per cell, $R^*$', fontsize=8)
_ = ax[0].set_ylabel(r'$\frac{\partial \Delta F }{ \partial R^*}$ [$k_BT$]',
    fontsize=8)
_ = ax[1].set_xlabel(r'DNA binding energy, $\Delta\varepsilon_{RA}$ [$k_BT$]',
    fontsize=8)
_ = ax[1].set_ylabel(r'$\frac{\partial \Delta F }{ \partial \Delta\varepsilon_{RA}^*}$ [$k_BT$]',
    fontsize=8)

# Special formatting. 
_ = ax[1].set_ylim([-1.5, -0.5])
for a in ax:
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)




plt.tight_layout()