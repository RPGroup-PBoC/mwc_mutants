# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
colors = mut.viz.personal_style()
constants = mut.thermo.load_constants()

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

# Define the reference state
R_0 = 260
epRA_0 = -13.9
c_0 = 50
c_smooth = np.logspace(-2, 4, 500)
c_data = c_smooth[::40]

# Compute the ground truth
ref_params = {'R':R_0, 'effector_conc':c_0,
             'ep_RA':epRA_0, 'Ka':constants['Ka'],
              'Ki': constants['Ki'], 
              'ep_AI':constants['ep_AI'],
              'n_sites':constants['n_sites']}

per_params = {'R':R_0, 'effector_conc':c_smooth,
             'ep_RA':epRA_0, 'Ka':constants['Ka'],
             'Ki':constants['Ki'], 'ep_AI':constants['ep_AI'],
             'n_sites':constants['n_sites']}

ground_truth = dbohr(ref_params, per_params)

ref = mut.thermo.SimpleRepression(R=R_0, effector_conc=c_0,
                                      ep_r=epRA_0, ka=constants['Ka'], 
                                      ki=constants['Ki'], n_sites=constants['n_sites'],
                                      ep_ai=constants['ep_AI'])

ref_arch = mut.thermo.SimpleRepression(R=R_0, effector_conc=c_smooth,
                                      ep_r=epRA_0, ka=constants['Ka'], 
                                      ki=constants['Ki'], n_sites=constants['n_sites'],
                                      ep_ai=constants['ep_AI'])
ref_bohr = -ref.bohr_parameter()

# Generate three data sets. 

# ####################################
# VARYING ∆Ep_RA
# ###################################
epRA_actual = [-17, -10]
ep, c = np.meshgrid(epRA_actual, c_data)
Ka_actual = constants['Ka']
Ki_actual = constants['Ki']
epRA_arch = mut.thermo.SimpleRepression(R=R_0, effector_conc=c,
                                      ep_r=ep, ka=constants['Ka'], 
                                      ki=constants['Ki'], n_sites=constants['n_sites'],
                                      ep_ai=constants['ep_AI'])
epRA_fc = epRA_arch.fold_change()
epRA_emp_df = ref_bohr - np.log((1/epRA_fc) -1)

# ###################################
# VARYING ∆EP_AI
# ###################################
epAI_actual =  [12, -2]
ep, c = np.meshgrid(epAI_actual, c_data)
Ka_actual = constants['Ka']
Ki_actual = constants['Ki']
epAI_arch = mut.thermo.SimpleRepression(R=R_0, effector_conc=c,
                                      ep_r=epRA_0, ka=constants['Ka'], 
                                      ki=constants['Ki'], n_sites=constants['n_sites'],
                                      ep_ai=ep)
epAI_fc = epAI_arch.fold_change()
epAI_emp_df = ref_bohr - np.log((1/epAI_fc) -1)


# ###################################
# VARYING Ka/Ki
# ###################################
KaKi_actual = np.array([1, 1E4]) 
Ka_actual = KaKi_actual * constants['Ki']
Ki_actual = constants['Ki']
ka, c = np.meshgrid(Ka_actual, c_data)
kaki_arch = mut.thermo.SimpleRepression(R=R_0, effector_conc=c,
                                      ep_r=epRA_0, ka=ka, 
                                      ki=constants['Ki'], n_sites=constants['n_sites'],
                                      ep_ai=constants['ep_AI'])
kaki_fc = kaki_arch.fold_change()
kaki_emp_df = ref_bohr - np.log((1/kaki_fc) -1)


# Instantiate the figure axis. 
fig, ax = plt.subplots(3, 3, figsize=(7, 6))

# Plot the reference induction profile. 
for i in range(3):
    ax[i, 0].plot(c_smooth / c_0, ref_arch.fold_change(), 'k-', lw=1, label='ref.')
    ax[i, 1].plot(c_smooth / c_0, ground_truth, '-', color='rebeccapurple', lw=1)
    ax[i, 2].plot([-10, 10], [-10, 10], ':', lw=1, alpha=0.5, color='k')

# ###################################
# Ep_RA 
# ###################################
# Plot the fake epRA data sets. 
ax[0, 0].plot(c_data / c_0, epRA_fc[:, 0], '.', color='firebrick', label=epRA_actual[0])
ax[0, 0].plot(c_data / c_0, epRA_fc[:, 1], '.', color=colors[3], label=epRA_actual[1])

# Plot the df v c/co
ax[0, 1].plot(c_data / c_0, epRA_emp_df[:, 0], '.', color='firebrick')
ax[0, 1].plot(c_data / c_0, epRA_emp_df[:, 1], '.', color=colors[3])

# Plot the self agreement. 
ax[0, 2].plot(ground_truth[::40], epRA_emp_df[:, 0], '.', color='firebrick')
ax[0, 2].plot(ground_truth[::40], epRA_emp_df[:, 1], '.', color=colors[3])
    
# ###################################
# Ep_AI
# ###################################
# Plot the fake eAI data sets. 
ax[1, 0].plot(c_data / c_0, epAI_fc[:, 0], '.', color='firebrick', label=epAI_actual[0])
ax[1, 0].plot(c_data / c_0, epAI_fc[:, 1], '.', color=colors[3], label=epAI_actual[1])

# Plot the df v c/co
ax[1, 1].plot(c_data / c_0, epAI_emp_df[:, 0], '.', color='firebrick')
ax[1, 1].plot(c_data / c_0, epAI_emp_df[:, 1], '.', color=colors[3])

# Plot the self agreement. 
ax[1, 2].plot(ground_truth[::40], epAI_emp_df[:, 0], '.', color='firebrick')
ax[1, 2].plot(ground_truth[::40], epAI_emp_df[:, 1], '.', color=colors[3])
    
# ###################################
# KaKi 
# ###################################
# Plot the fake eAI data sets. 
ax[2, 0].plot(c_data / c_0, kaki_fc[:, 0], '.', color='firebrick', label='1')
ax[2, 0].plot(c_data / c_0, kaki_fc[:, 1], '.', color=colors[3], label='$10^4$')

# Plot the df v c/co
ax[2, 1].plot(c_data / c_0, kaki_emp_df[:, 0], '.', color='firebrick')
ax[2, 1].plot(c_data / c_0, kaki_emp_df[:, 1], '.', color=colors[3])

# Plot the self agreement. 
ax[2, 2].plot(ground_truth[::40], kaki_emp_df[:, 0], '.', color='firebrick')
ax[2, 2].plot(ground_truth[::40], kaki_emp_df[:, 1], '.', color=colors[3])
     
# Formating and labeling
row_labels = [r'varying $\Delta\varepsilon_{RA}$',
             r'varying $\Delta\varepsilon_{AI}$',
             'varying $K_A / K_I$']
leg_titles = [r'$\Delta\varepsilon_{RA}$ [$k_BT$]', r'$\Delta\varepsilon_{AI}$ [$k_BT$]', '$K_A / K_I$']
for i in range(3):
    ax[i,0].text(-.46, 0.72, row_labels[i], rotation='vertical', 
                 backgroundcolor='#f1f2f6', transform=ax[i,0].transAxes, fontsize=8)
    ax[i, 0].set_ylim([0, 1.4])
    ax[i, 2].set_ylim([-10, 10])
    ax[i, 2].set_xlim([-10, 10])
    ax[i, 1].set_ylim([-10, 10])
    ax[i, 0].set_xscale('log')
    ax[i, 0].set_ylabel('fold-change')
    ax[i, 0].set_xlabel('$c / c_0$')
    ax[i, 1].set_xlabel('$c / c_0$')
    ax[i, 1].set_xscale('log')
    ax[i, 1].set_ylabel(r'$\Delta F$ [$k_BT$]')
    ax[i, 2].set_ylabel(r'measured $\Delta F$ [$k_BT$]')
    ax[i, 2].set_xlabel(r'predicted $\Delta F$ [$k_BT$]')
    
    # Set the labels
    _leg = ax[i, 0].legend(fontsize=8, handlelength=1, title=leg_titles[i], loc='upper left')
    _leg.get_title().set_fontsize(8)
    
# Add panel labels. 
fig.text(0.02, 0.99, '(a)', fontsize=8)
fig.text(0.02, 0.66, '(b)', fontsize=8)
fig.text(0.02, 0.33, '(c)', fontsize=8)
plt.tight_layout()


plt.savefig('./example_deviations.pdf')
