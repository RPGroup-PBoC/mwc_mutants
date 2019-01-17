# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mut.thermo
import mut.viz
constants = mut.thermo.load_constants()
colors = mut.viz.personal_style()
# Set the ranges. 
c_range = np.logspace(-2, 4, 200)
R_range = np.logspace(0, 4, 200)
epRA_range = np.linspace(-20, -5, 200)

# Define ec50 (approx)
c_0 = 50 # in µM
R_0 = 260
epRA_0 = -13.9 # in kBT

# Define the delta bohr function
def dbohr(ref_params, per_params):
    # Compute the pact for reference and perturbed
    ref_pact = mut.thermo.MWC(effector_conc=ref_params['effector_conc'],
                             ka=ref_params['Ka'], ki=ref_params['Ki'],
                             ep_ai=ref_params['ep_AI'], 
                             n_sites=ref_params['n_sites']).pact()
    per_pact = mut.thermo.MWC(effector_conc=per_params['effector_conc'],
                             ka=per_params['Ka'], ki=per_params['Ki'],
                             ep_ai=per_params['ep_AI'], 
                             n_sites=per_params['n_sites']).pact()
    
    # Compute and return the delta Bohr. 
    delta_bohr = np.log(ref_pact/per_pact) +\
                 np.log(ref_params['R']/per_params['R']) -\
                 (ref_params['ep_RA'] - per_params['ep_RA'])
    return delta_bohr


# #####################################
# COMPUTE SURFACES
# ####################################
cr_mesh, rc_mesh = np.meshgrid(c_range, R_range)
cep_mesh, epc_mesh = np.meshgrid(c_range, epRA_range)
epr_mesh, rep_mesh = np.meshgrid(epRA_range, R_range)
ref_params = {'effector_conc': c_0,
             'Ka':constants['Ka'],
             'Ki':constants['Ki'],
             'ep_AI': constants['ep_AI'],
             'n_sites': constants['n_sites'],
             'ep_RA': epRA_0, 'R':R_0}
per_c_epRA_params = {'effector_conc': cep_mesh,
                  'Ka':constants['Ka'],
                  'Ki':constants['Ki'],
                  'ep_AI':constants['ep_AI'],
                  'n_sites':constants['n_sites'],
                  'ep_RA':epc_mesh, 'R':R_0}
per_c_R_params = {'effector_conc': cr_mesh,
                  'Ka':constants['Ka'],
                  'Ki':constants['Ki'],
                  'ep_AI':constants['ep_AI'],
                  'n_sites':constants['n_sites'],
                  'ep_RA':epRA_0, 'R':rc_mesh}
per_ep_R_params = {'effector_conc': c_0,
                  'Ka':constants['Ka'],
                  'Ki':constants['Ki'],
                  'ep_AI':constants['ep_AI'],
                  'n_sites':constants['n_sites'],
                  'ep_RA':epr_mesh, 'R':rep_mesh}
c0_r = dbohr(ref_params, per_c_R_params)
c0_ep = dbohr(ref_params, per_c_epRA_params)
r_ep = dbohr(ref_params, per_ep_R_params)


# ####################################
# PLOTTING
# ####################################
fig, ax = plt.subplots(1, 3, figsize=(7, 3))
for a in ax:
    a.grid(False)
   
# Plot the heat maps
ax[0].imshow(c0_r, vmin=-10, vmax=12, 
            interpolation='none', cmap='bone')
ax[1].imshow(c0_ep, vmin=-10, vmax=12, 
            interpolation='none', cmap='bone')
ax[2].imshow(r_ep, vmin=-10, vmax=12, 
            interpolation='none', cmap='bone')

# Plot the contours. 
c1 = ax[0].contour(c0_r, levels=[-5, 0, 5], linestyles='-', colors='w',
                  linewidths=1)
c2 = ax[1].contour(c0_ep, levels=[-5, 0, 5], linestyles='-', colors='w',
                  linewidths=1)
c3 = ax[2].contour(r_ep, levels=[-5, 0, 5], linestyles='-', colors='w', 
                  linewidths=1)

# Add contour labels. 
cl0 = ax[0].clabel(c1, c1.levels, inline=True, fmt='$\,$ %d $k_BT$', fontsize=8,
            manual=[(100, 150), (100, 100), (150, 25)], inline_spacing=3)
cl1 = ax[1].clabel(c2, c2.levels, inline=True, fmt='$\,$ %d $k_BT$', fontsize=8,
            manual=[(50, 75), (50, 125), (100, 150)], inline_spacing=3)
cl2 = ax[2].clabel(c3, c3.levels, inline=True, fmt='$\,$ %d $k_BT$', fontsize=8,
            manual=[(20, 150), (100, 100), (150, 50)], inline_spacing=3)
for l in cl0 + cl1 + cl2:
    l.set_rotation(0)
for a in ax:
    a.spines['top'].set_visible(True)
    a.spines['right'].set_visible(True)
    
# Format the labeling
ep_ticks = epRA_range[::40]
ax[0].set_xticks(np.arange(0, 200, 60))
ax[0].set_xticklabels(['$10^{-2}$', '$10^{-1}$',
                      '$10^0$', '$10^1$'])
ep_ticks = list(np.round(ep_ticks))
ep_ticks.reverse()
ep_ticks.insert(0, '')
ax[0].set_yticklabels(ep_ticks)

ax[1].set_xticks(np.arange(0, 200, 60))
ax[1].set_xticklabels(['$10^{-2}$', '$10^{-1}$',
                      '$10^0$', '$10^1$'])
ax[1].set_yticklabels(['', '$10^4$', '$10^3$', '$10^2$', '$10^1$', '$10^0$']) 
ax[2].set_xticklabels(['', '$10^0$', '$10^1$', '$10^2$', '$10^3$', '$10^4$']) 
ax[2].set_yticklabels(ep_ticks)

# Add axis labels. 
ax[0].set_ylabel(r'$\Delta\varepsilon_{RA}$ [$k_BT$]')
ax[0].set_xlabel(r'$(c / c_0)$')
ax[1].set_xlabel(r'$(c / c_0)$')
ax[1].set_ylabel('repressors per cell')
ax[2].set_xlabel('repressors per cell')
ax[2].set_ylabel(r'$\Delta\varepsilon_{RA}$ [$k_BT$]')

# Add panel labels. 
fig.text(0, 0.8, '(a)', fontsize=8)
fig.text(0.34, 0.8, '(b)', fontsize=8)
fig.text(0.66, 0.8, '(c)', fontsize=8)
plt.tight_layout()
plt.savefig('./deltaF_heatmaps.pdf', bbox_inches='tight')



