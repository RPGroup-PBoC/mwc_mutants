# -*- coding: utf-8 -*- import sys sys.path.insert(0, '../../')
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import mut.thermo
import mut.stats
import mut.viz
pboc = mut.viz.color_selector('pboc')
constants = mut.thermo.load_constants()
mut.viz.plotting_style()

# Define the reference states
r0 = 200 # per cell
epRA0 = -15 # in kBT
ka0 = 200
ki0 = 0.5
epAI0 = 5
c_range = np.logspace(-2, 4, 500)

ref_pact = mut.thermo.MWC(effector_conc=c_range,
                         ka=ka0, ki=ki0, ep_ai=epAI0,
                         n_sites=constants['n_sites']).pact()
ref_stats = {'pact':ref_pact, 'R':r0, 'ep_RA':epRA0}

# Compute the changed states. 
r_range = [1, 10, 100, 1000]
epRA_range = [-25, -20, -10, -5]
ka_range = [1, 100, 500, 1000]
ki_range = [0.5, 0.5, 0.5, 0.5]
epAI_range = [-2, 0, 1, 10]

# Define a function to compute the delta bohr
def delta_bohr(ref_stats, mut_stats):
    pact_factor = -np.log(ref_stats['pact']/ mut_stats['pact'])
    rep_factor = -np.log(ref_stats['R'] / mut_stats['R'])
    ep_factor = ref_stats['ep_RA'] - mut_stats['ep_RA']
    return pact_factor + rep_factor + ep_factor

r, c = np.meshgrid(r_range, c_range)
ref_pact = mut.thermo.MWC(effector_conc=c,
                         ka=ka0, ki=ki0, ep_ai=epAI0,
                         n_sites=constants['n_sites']).pact()
ref_stats = {'pact':ref_pact, 'R':r0, 'ep_RA':epRA0}
r_c_stats = {'pact':ref_pact, 'R':r, 'ep_RA':epRA0}
r_c_dbohr = delta_bohr(ref_stats, r_c_stats)

# Set up the figure canvas
fig, ax = plt.subplots(2, 4, figsize=(5.5, 2.5))

for i in range(len(r_range)):
    ax[1, 0].plot(c, r_c_dbohr[:, i])
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_xscale('symlog')
    a.hlines(0, -0.1, 1E4, pboc['dark_green'], linestyle='--')
    a.set_ylim([-5, 5])
# plt.tight_layout()


# ref_stats = {'R':r0, 'pact':pact0, 'ep_RA':epRA0}

# # Define the pact
# pact0 = 0.5
    
# # Define the parameter ranges.
# r_range = np.logspace(0, 3, 500)
# c_range = np.logspace(-2, 4, 500)
# pact_range = mut.thermo.MWC(ka=constants['Ka'], ki=constants['Ki'],
#                       ep_ai=constants['ep_AI'], effector_conc=c_range,
#                       n_sites=constants['n_sites']).pact()
# epRA_range = np.linspace(-20, -5, 500)

# # Define the delta F function
# def delta_bohr(ref_stats, mut_stats):
#     pact_factor = -np.log(ref_stats['pact']/ mut_stats['pact'])
#     rep_factor = -np.log(ref_stats['R'] / mut_stats['R'])
#     ep_factor = ref_stats['ep_RA'] - mut_stats['ep_RA']
#     return pact_factor + rep_factor + ep_factor


# # Mesh the parameter pairs and compute the delta F

# # Repressors v epRA
# r_ep, ep_r = np.meshgrid(r_range, epRA_range)
# r_ep_stats = {'R':r_ep, 'pact':pact0, 'ep_RA':ep}
# r_ep_dbohr = delta_bohr(ref_stats, r_ep_stats)

# # Repressors v pact 
# r_p, p_r = np.meshgrid(r_range, pact_range)
# r_p_stats = {'R':r_p, 'pact':p_r, 'ep_RA':epRA0}
# r_p_dbohr = delta_bohr(ref_stats, r_p_stats)

# # epRA v pact 
# ep_p, p_ep = np.meshgrid(epRA_range, pact_range)
# ep_p_stats = {'R':r0, 'pact':p_ep, 'ep_RA':ep_p}
# ep_p_dbohr = delta_bohr(ref_stats, ep_p_stats)


# fig, ax = plt.subplots(1, 3, figsize=(5.5, 2.5))

# ax[0].imshow(r_ep_dbohr, cmap='bone', origin='lower', vmin=-14, vmax=16)
# cs0 = ax[0].contour(r_ep_dbohr, levels=[-10,-5, 0, 5, 10], colors='w',
#              linestyles='-', origin='lower')
# ax[0].clabel(cs0, colors=['w'], inline=True, fmt='%1.0f $k_BT$', rightside_up=True,
#             fontsize=8)

# ax[1].imshow(r_p_dbohr, cmap='bone', origin='lower', vmin=-14, vmax=16)
# cs1 = ax[1].contour(r_p_dbohr, levels=[-8, -5, -2.5, 0, 2.5, 5], colors='w',
#              linestyles='-', origin='lower')
# ax[1].clabel(cs1, colors=['w'], inline=True, fmt='%1.0f $k_BT$', rightside_up=True,
#             fontsize=8)

# ax[2].imshow(ep_p_dbohr, cmap='bone', origin='lower', vmin=-14, vmax=16)
# cs2 = ax[2].contour(ep_p_dbohr, levels=[-10, -5, 0, 3], colors='w',
#              linestyles='-', origin='lower')
# ax[2].clabel(cs2, colors=['w'], inline=True, fmt='%1.0f $k_BT$', rightside_up=True,
#             fontsize=8)

# for a in ax:
#     a.grid(False)
#     a.set_frame_on(True)
#     a.xaxis.set_tick_params(labelsize=8)
#     a.yaxis.set_tick_params(labelsize=8)
# plt.tight_layout()    
    