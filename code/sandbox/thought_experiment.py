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
c_smooth[0] = 0
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
epRA_actual = np.array([-17, -10])
ep, c = np.meshgrid(epRA_actual, c_smooth)
Ka_actual = constants['Ka']
Ki_actual = constants['Ki']
epRA_arch = mut.thermo.SimpleRepression(R=R_0, effector_conc=c,
                                      ep_r=ep, ka=constants['Ka'], 
                                      ki=constants['Ki'], n_sites=constants['n_sites'],
                                      ep_ai=constants['ep_AI'])
epRA_fc = epRA_arch.fold_change()
epRA_emp_df = ref_bohr - np.log((1/epRA_fc) -1)

# compute the theoretical. ddf. 
epRA_ddf = np.array([ground_truth - epRA_emp_df[:, 0],
                     ground_truth - epRA_emp_df[:, 1]]).T

# ###################################
# VARYING ∆EP_AI
# ###################################
epAI_actual =  [8, -2]
ep, c = np.meshgrid(epAI_actual, c_smooth)
Ka_actual = constants['Ka']
Ki_actual = constants['Ki']
epAI_arch = mut.thermo.SimpleRepression(R=R_0, effector_conc=c,
                                      ep_r=epRA_0, ka=constants['Ka'], 
                                      ki=constants['Ki'], n_sites=constants['n_sites'],
                                      ep_ai=ep)
epAI_fc = epAI_arch.fold_change()
epAI_emp_df = ref_bohr - np.log((1/epAI_fc) -1)
epAI_ddf = np.array([ground_truth - epAI_emp_df[:, 0],
                    ground_truth - epAI_emp_df[:, 1]]).T

# ###################################
# VARYING Ka/Ki
# ###################################
KaKi_actual = np.array([0.1,  1000 ]) 
Ka_actual = KaKi_actual * constants['Ki']
Ki_actual = constants['Ki'] 
ka, c = np.meshgrid(Ka_actual, c_smooth)
kaki_arch = mut.thermo.SimpleRepression(R=R_0, effector_conc=c,
                                      ep_r=epRA_0, ka=ka, 
                                      ki=Ki_actual, n_sites=constants['n_sites'],
                                      ep_ai=constants['ep_AI'])
kaki_fc = kaki_arch.fold_change()
kaki_emp_df = ref_bohr - np.log((1/kaki_fc) -1)
kaki_ddf = np.array([ground_truth - kaki_emp_df[:, 0],
                    ground_truth - kaki_emp_df[:, 1]]).T

# #####################################
# PLOTTING 
# ######################################
def plotting_fn(c, c_0, fc, emp_df, 
                ddf, ref_arch=ref_arch, ground_truth=ground_truth,  
                skip=38, labels=['__nolegend__', '__nolegend__'], 
                title=None, savename=None):

    fig, ax = plt.subplots(1, 3, figsize=(7, 2), dpi=100)

    # Plot the reference induction profile. 
    ax[0].plot(c / c_0, ref_arch.fold_change(), 'k-', lw=1, label='__nolegend__')
    ax[1].plot(c / c_0, ground_truth, 'k-', lw=1, label='__nolegend__')
    ax[2].plot([-1E3, 1E3], [0, 0], ':', lw=1, alpha=0.5, color='k', label='__nolegend__')
        
    # Plot the two generated data sets. 
    ax[0].plot(c[::skip]/ c_0, fc[::skip, 0], '.', color='dodgerblue', ms=4, 
              label=labels[0])
    ax[0].plot(c / c_0, fc[:, 0], '-', 
               color='dodgerblue', ms=4, lw=0.5, alpha=0.75,  label='__nolegend__')
    ax[0].plot(c[::40] / c_0, fc[::40, 1], '.', color='firebrick', ms=4,
              label=labels[1])
    ax[0].plot(c / c_0, fc[:, 1], '-', 
               color='firebrick', ms=4, lw=0.5, alpha=0.75, label='__nolegend__')
    
    # Plot the anticipated collapse. 
    ax[1].plot(c[::skip] / c_0, emp_df[::skip, 0], '.', color='dodgerblue', ms=4,
              label=labels[0])
    ax[1].plot(c / c_0, emp_df[:,  0], '-',
              color='dodgerblue', lw=0.5, alpha=0.75,  label='__nolegend__')
    ax[1].plot(c[::skip] / c_0, emp_df[::skip, 1], '.', color='firebrick', ms=4,
              label=labels[1])
    ax[1].plot(c / c_0, emp_df[:,  1], '-',
              color='firebrick', lw=0.5, alpha=0.75, label='__nolegend__')
    
    # Plot the ddf
    ax[2].plot(c[::skip] / c_0, ddf[::skip, 0], '.', color='dodgerblue', ms=4,
              label=labels[0])
    ax[2].plot(c / c_0, ddf[:, 0], '-', lw=0.5, color='dodgerblue', alpha=0.75,
              label='__nolegend__')
    ax[2].plot(c[::skip] / c_0, ddf[::skip, 1], '.', color='firebrick', ms=4,
              label=labels[1])
    ax[2].plot(c / c_0, ddf[:, 1], '-', lw=0.5, color='firebrick', alpha=0.75,
              label='__nolegend__')
    
    # Format the axes. 
    for i in range(3):
        ax[i].set_xscale('symlog', linthreshx=c_smooth[1])
    
    # Add labels. 
    ax[0].set_xlabel('$c / c_0$')
    ax[0].set_ylabel('fold-change')
    ax[0].set_ylim([-0.2, 1.2])
    ax[1].set_xlabel('$c / c_0$')
    ax[1].set_ylabel(r'$\Delta F$ [$k_BT$]')
    ax[1].set_ylim([-10,10])
    ax[2].set_xlim([0, 2E2])
    ax[2].set_ylim([-10, 10])
    ax[2].set_xlabel(r'$c / c_0$')
    ax[2].set_ylabel(r'$\Delta F^* - \Delta F$ [$k_BT$]')
    
    handles, _ = ax[0].get_legend_handles_labels()
    ax[0].text(-0.45, 1.1, '(a)', fontsize=8, transform=ax[0].transAxes)
    ax[1].text(-0.4, 1.1, '(b)', fontsize=8, transform=ax[1].transAxes)
    ax[2].text(-0.4, 1.1, '(c)', fontsize=8, transform=ax[2].transAxes)
    plt.subplots_adjust(wspace=0.5)
    if savename != None:
        plt.savefig(savename, bbox_inches='tight')
    return [fig, ax, leg]


_ = plotting_fn(c_smooth, c_0, epRA_fc, epRA_emp_df, epRA_ddf, 
                                savename='epRA_ddf.pdf')
_ = plotting_fn(c_smooth, c_0, kaki_fc, kaki_emp_df, kaki_ddf, 
                                savename='kaki_ddf.pdf')


