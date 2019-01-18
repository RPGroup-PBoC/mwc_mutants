# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
import imp
imp.reload(mut.viz)
colors = mut.viz.personal_style()
constants = mut.thermo.load_constants()

# Define the deltaF of perfect agreement
deltaF_range = np.linspace(-8, 8, 200)
n_meas = 10 

# Compute the various couplings
deltaF_meas = np.linspace(-4, 4, n_meas)
deltaF_add = deltaF_meas + 2
deltaF_mult = deltaF_meas * 0.5
deltaF_comp = deltaF_meas**2  / 5

# Define the range of delta 
fig, ax = plt.subplots(2, 2, figsize=(5, 4.5), sharex=True, sharey=True)
ax = ax.ravel()
for a in ax:
    a.plot(deltaF_range, deltaF_range, '--', color='rebeccapurple', lw=0.5)
ax[0].plot(deltaF_meas, deltaF_meas, '.', color=colors[3], ms=5)
ax[1].plot(deltaF_meas, deltaF_add, '.', color=colors[3], ms=5)
ax[2].plot(deltaF_meas, deltaF_mult, '.', color=colors[3], ms=5)
ax[3].plot(deltaF_meas, deltaF_comp, '.', color=colors[3], ms=5)
