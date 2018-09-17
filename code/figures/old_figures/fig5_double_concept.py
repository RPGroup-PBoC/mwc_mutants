# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.thermo
import mut.viz
mut.viz.plotting_style()
colors = mut.viz.color_selector('pboc')


# Set the constants.
ep_ai = 4.5
ep_r_wt = -12
ka_wt = 150E-6
ki_wt = 0.5E-6
ka_mut = 100E-6
ki_mut = .1E-6
ep_r_mut = -18
R = 300

# Mesh the parameter values 
c_range =  np.logspace(-8, -2, 500)
ep_r = [ep_r_wt, ep_r_mut]
_ka = [ka_wt, ka_mut]
_ki = [ki_wt, ki_mut]
c, ep, ka, ki = np.meshgrid(c_range, ep_r, _ka, _ki)

# Set up the architecture.
arch = mut.thermo.SimpleRepression(R, ep, effector_conc=c, ka=ka, ki=ki, ep_ai=ep_ai)
fc = arch.fold_change()

# Define the colors 
color_key = [colors['blue'], colors['red'], colors['purple']]

# %% Set up the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(3, 3)) 
# ax[0].axis('off')
# ax = ax[1]

# Set up the axis labels
ax.set_xlabel('$c / K_A^{wt}$', fontsize=8)
ax.set_ylabel('fold-change', fontsize=8)
ax.set_xscale('log')
ax.tick_params(labelsize=8)
ax.set_ylim([-0.1, 1.05])
ax.set_xlim([1E-8/ka_wt, 1E-2/ka_wt])
ax.grid(False)

# Plot the fold change.
_ = ax.plot(c_range / ka_wt, fc[0, :, 1, 1], label='inducer\nmutant',
    color = color_key[0])
_ = ax.plot(c_range / ka_wt, fc[1, :, 0, 0],  label='DNA\nmutant',
    color = color_key[1])
_ = ax.plot(c_range / ka_wt, fc[1, :, 1, 1], label='double\nmutant',
    color = color_key[2])
_ = ax.plot(c_range / ka_wt, fc[0, :, 0, 0], label='wild-type',
    color='slategray', alpha=1)
_ = ax.legend(fontsize=8, handlelength=1)
plt.tight_layout()
plt.savefig('fig5_plots.svg', bbox_inches='tight')


