# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../../')
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import mut.viz
import mut.thermo
constants = mut.thermo.load_constants()
colors = mut.viz.personal_style()


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

# Define the reference state. 
epRA_0 = -14
Ka_0 = 200
Ki_0 = 0.5
epAI_0 = 5
R_0 = 200
n_sites = 2
c_0 = np.logspace(-2, 4,  500)
ref_arch = mut.thermo.SimpleRepression(R=R_0, effector_conc=c_0, ka=Ka_0, ki=Ki_0, 
                                   ep_ai=epAI_0, ep_r=epRA_0, n_sites=n_sites) 
ref_fc = ref_arch.fold_change()
ref_params = {'effector_conc':c_0, 'Ka':Ka_0,  'Ki':Ki_0, 'ep_RA':epRA_0,
              'ep_AI':epAI_0, 'R':R_0, 'n_sites':n_sites}
delF = np.zeros(len(c_0))

# Define the perturbed states. 
epRA = [-17, -10]
c, ep = np.meshgrid(c_0, epRA)
mut_arch = mut.thermo.SimpleRepression(R=R_0, effector_conc=c, ka=Ka_0, ki=Ki_0, 
                                   ep_ai=epAI_0, ep_r=ep, n_sites=n_sites) 
mut_fc = mut_arch.fold_change()
mut_params = {'effector_conc':c_0, 'Ka':Ka_0, 'Ki':Ki_0, 'ep_RA':epRA[0],
               'ep_AI':epAI_0, 'R':R_0, 'n_sites':n_sites}

mut_delF = np.zeros((2, len(c_0)))
for i, ep in enumerate(epRA):
    mut_params['ep_RA'] = ep
    mut_delF[i, :] = dbohr(ref_params, mut_params)

# Set up the figure. 
fig, ax = plt.subplots(1, 2, figsize=(6, 3))

# Plot the fold-change of the reference state
ax[0].plot(c_0/ Ka_0, ref_fc, '-', color=colors[0],  label='reference')
ax[0].plot(c_0[::50] / Ka_0, ref_fc[::50], 'o', color=colors[0], label='__nolegend__')

# Plot the mutant fold-changes
ax[0].plot(c_0 / Ka_0, mut_fc[0, :], '-', color=colors[2],  label='mutant 1')
ax[0].plot(c_0[::50]/ Ka_0, mut_fc[0, ::50], 'o', color=colors[2], label='__nolegend__')
ax[0].plot(c_0/ Ka_0, mut_fc[1, :], '-', color=colors[4],  label='mutant 2')
ax[0].plot(c_0[::50] / Ka_0, mut_fc[1, ::50], 'o', color=colors[4], label='__nolegend__')

# Plot the delta Bohr. 
ax[1].plot(c_0 / Ka_0, delF, '-', color=colors[0], label='reference')
ax[1].plot(c_0[::50] / Ka_0, delF[::50], 'o', color=colors[0], label='__nolegend__')
ax[1].plot(c_0 / Ka_0, mut_delF[0, :], '-', color=colors[2], label='mutant 1')
ax[1].plot(c_0[::50] / Ka_0, mut_delF[0, ::50], 'o', color=colors[2], label='__nolegend__')
ax[1].plot(c_0 / Ka_0, mut_delF[1, :], '-', color=colors[4], label='mutant 2')
ax[1].plot(c_0[::50] / Ka_0, mut_delF[1, ::50], 'o', color=colors[4], label='__nolegend__')

# Format the axes
for a in ax:
    a.set_xscale('log')
ax[0].set_ylim([-0.15, 1.15])
ax[1].set_ylim([-5, 5])

# Label the axes
ax[0].set_xlabel('$c / K_{A_0}$')
ax[1].set_xlabel('$c / K_{A_0}$')
ax[0].set_ylabel('fold-change')
ax[1].set_ylabel('$\Delta F$ [$k_BT$]')

# Add a legend and panel labels
ax[0].legend(handlelength=1)
fig.text(0.0, 1, '(a)', fontsize=9)
fig.text(0.5, 1, '(b)', fontsize=9)
plt.tight_layout()
plt.savefig('thought_experiment_DNA_mutant.pdf', bbox_inches='tight')