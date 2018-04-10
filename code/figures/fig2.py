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
mut_ka = 30E-6  # in M
mut_ki = 3E-6  # in M
c_range = np.logspace(-8, -1, 300)
ep_range = np.linspace(-20, -5, 300)
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

# Example props.
leak = arch.leakiness()
dyn_rng = arch.dynamic_range()


# %%
# Set up the figure canvas.
fig = plt.figure(figsize=(6, 6))
gs = gridspec.GridSpec(8, 10)

# Part A
ax0 = fig.add_subplot(gs[0:4, 0:1])
ax1 = fig.add_subplot(gs[0:4, 1:6])
ax2 = fig.add_subplot(gs[0:2, 6:])
ax3 = fig.add_subplot(gs[2:4:, 6:])

# Part B
ax4 = fig.add_subplot(gs[4:, 0:1])
ax5 = fig.add_subplot(gs[4:, 1:6])
ax6 = fig.add_subplot(gs[4:6, 6:])
ax7 = fig.add_subplot(gs[6::, 6:])

# Format the axes
ax = [ax0, ax1, ax2, ax3, ax4, ax5, ax6, ax7]
for a in ax:
    a.set_yscale('log')

ka_ki = c_range / wt_ka
for a in [ax1, ax5]:
    a.set_xlim([ka_ki[0], ka_ki[-1]])

for a in [ax2, ax3]:
    a.set_xlim([-20, -5])
for a in [ax6, ax7]:
    a.set_xlim([(ka_range / wt_ki)[0], (ka_range / wt_ki)[-1]])

ax[6].set_xscale('log')
ax[7].set_xscale('log')
ax[1].set_xscale('log')
ax[5].set_xscale('log')
ax[0].axis('off')
ax[4].axis('off')

# Add labels.
ax[1].set_xlabel('$c / K_A^{wt}$')
ax[5].set_xlabel('$c / K_A^{wt}$')
ax[1].set_ylabel('fold-change')
ax[5].set_ylabel('fold-change')
ax[2].set_xlabel('DNA binding energy [$k_BT$]')
ax[3].set_xlabel('DNA binding energy [$k_BT$]')
ax[2].set_ylabel('leakiness')
ax[3].set_ylabel('dynamic\n range')
ax[6].set_xlabel('$K_A / K_I$')
ax[7].set_xlabel('$K_A / K_I$')
ax[6].set_ylabel('leakiness')
ax[7].set_ylabel('dynamic\n range')

# PLot the titration curves.
_ = ax1.plot(c_range / ka[0], fc[0, 0, 0, :], '-', color=colors['red'],
             lw=1.5, label='wild-type')
_ = ax1.plot(c_range / ka[0], fc[0, 1, 0, :], '-', color=colors['blue'],
             lw=1.5, label='mutant')

_ = ax5.plot(c_range / ka[0], fc[0, 0, 0, :], '-', color=colors['red'],
             lw=1.5, label='wild-type')
_ = ax5.plot(c_range / ka[0], fc[1, 0, 1, :], '-', color=colors['green'],
             lw=1.5, label='mutant')

ax1.legend(handlelength=1)
ax5.legend(handlelength=1)
# Plot the properties from theory and examples..
_ = ax2.plot(ep_range, DNA_leak, 'k-', lw=1)
_ = ax3.plot(ep_range, DNA_dyn_rng, 'k-', lw=1)
_ = ax6.plot(ka_range/wt_ki, np.ones_like(ka_range) * IND_leak, 'k-', lw=1)
_ = ax7.plot(ka_range/wt_ki, IND_dyn_rng, 'k-', lw=1)

# Plot points from the example mutants.
_ = ax2.plot(wt_epR, leak[0, 0, 0, 0], 'o', color=colors['red'], ms=5)
_ = ax2.plot(mut_epR, leak[0, 1, 0, 0], 'o', color=colors['blue'], ms=5)
_ = ax3.plot(wt_epR, dyn_rng[0, 0, 0, 0], 'o', color=colors['red'], ms=5)
_ = ax3.plot(mut_epR, dyn_rng[0, 1, 0, 0], 'o', color=colors['blue'], ms=5)
_ = ax6.plot(wt_ka / wt_ki, leak[0, 0, 0, 0], 'o', color=colors['red'], ms=5)
_ = ax6.plot(mut_ka / mut_ki, leak[1, 0, 1, 0], 'o', color=colors['green'], ms=5)
_ = ax7.plot(wt_ka / wt_ki, dyn_rng[0, 0, 0, 0], 'o', color=colors['red'], ms=5)
_ = ax7.plot(mut_ka / mut_ki, dyn_rng[1, 0, 1, 0], 'o', color=colors['green'], ms=5)

# Format and save
plt.tight_layout(h_pad=0.6)

# Add panel labels.
fig.text(-0.3, 1.0, '(A)', fontsize=8)
fig.text(-0.3, 0.5, '(B)', fontsize=8)


plt.savefig('figures/fig{}_plots.pdf'.format(FIG_NO))

