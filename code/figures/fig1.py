# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
sys.path.insert(0, '../../')
import mut.viz
import mut.thermo
colors = mut.viz.pub_style()

FIG_NO = 2

# Define the architectural parameters
R = 300
wt_epR = -15  # in k_BT
wt_ka = 200E-6  # in M
wt_ki = 0.5E-6  # in M
mut_epR = -10  # in k_BT
mut_ka = 500E-6  # in M
mut_ki = 3E-6  # in M
c_range = np.logspace(-8, -1, 300)
ep_range = np.linspace(-20, -10, 300)
ka_range = np.logspace(-6, -3, 300)
ki_range = np.logspace(-3, 0, 300)

# Mesh together the parameters for evaluation.
ep_R = [wt_epR, mut_epR]
ka = [wt_ka, mut_ka]
ki = [wt_ki, mut_ki]

ep_mesh, ka_mesh, ki_mesh, c_mesh = np.meshgrid(ep_R, ka, ki, c_range)

# Instantiate the architecture
arch = mut.thermo.SimpleRepression(R, ep_r=ep_mesh, effector_conc=c_mesh, ka=ka_mesh,
                                   ki=ki_mesh, ep_ai=4.5, n_sites=2)

fc = arch.fold_change()

# Compute the properties for the dynamic range, satuation, and leakiness
ep_mesh, ka_mesh, ki_mesh = np.meshgrid(ep_range, ka_range, ki_range)


DNA_props = mut.thermo.SimpleRepression(R, ep_range, ka=wt_ka, ki=wt_ki, ep_ai=4.5, n_sites=2, effector_conc=0)
IND_props = mut.thermo.SimpleRepression(R, wt_epR, ka=ka_range, ki=wt_ki, ep_ai=4.5, n_sites=2, effector_conc=0)
DNA_leak = DNA_props.leakiness()
DNA_dyn_rng = DNA_props.dynamic_range()
IND_leak = IND_props.leakiness()
IND_dyn_rng = IND_props.dynamic_range()

# %%
# Set up the figure canvas.
fig = plt.figure(figsize=(6, 6))
gs = gridspec.GridSpec(8, 10)

# Part A
ax0 = fig.add_subplot(gs[0:4, 0:3])
ax1 = fig.add_subplot(gs[0:4, 2:6])
ax2 = fig.add_subplot(gs[0:2, 6:])
ax3 = fig.add_subplot(gs[2:4:, 6:])

# Part B
ax4 = fig.add_subplot(gs[4:, 0:3])
ax5 = fig.add_subplot(gs[4:, 2:6])
ax6 = fig.add_subplot(gs[4:6, 6:])
ax7 = fig.add_subplot(gs[6::, 6:])

axes = [ax0, ax1, ax2, ax3, ax4, ax5, ax6, ax7]

    # a.set_ylim([-0.1, 1.2])

ax2.set_yscale('log')
ax6.set_yscale('log')
ax7.set_yscale('log')
ax6.set_xscale('log')
ax7.set_xscale('log')
# Tie together the axes and format as necessary
ax0.axis('off')
ax4.axis('off')
ax1.set_xscale('log')
ax5.set_xscale('log')

# PLot the titration curves.
_ = ax1.plot(c_range / ka[0], fc[0, 0, 0, :], '-', color=colors['blue'],
             lw=2)
_ = ax1.plot(c_range / ka[0], fc[0, 1, 0, :], '-', color=colors['green'],
             lw=2)

_ = ax5.plot(c_range / ka[0], fc[0, 0, 0, :], '-', color=colors['blue'],
             lw=2)
_ = ax5.plot(c_range / ka[0], fc[1, 0, 1, :], '-', color=colors['green'],
             lw=2)


# Plot the properties.
_ = ax2.plot(ep_range, DNA_leak)
_ = ax3.plot(ep_range, DNA_dyn_rng) 

_ = ax6.plot(ka_range/wt_ki, np.ones_like(ka_range) * IND_leak)
_ = ax7.plot(ka_range/wt_ki, IND_dyn_rng)


plt.tight_layout()

