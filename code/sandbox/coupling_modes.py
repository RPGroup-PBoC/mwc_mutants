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
n_meas = 200

# Compute the various couplings
deltaF_meas = np.linspace(-8, 8, 200)
deltaF_linear = deltaF_meas + 5
deltaF_mult = 0.25 * deltaF_meas
deltaF_comp = 2 * np.cos(deltaF_meas)+1

# Generate a fake data set. 
c_range_data = np.logspace(-2, 4, 10)
c_range = np.logspace(-2, 4, 100)
ep_RA_pred = -15
ep_RA_actual = -13
ep_RA_ref = -15
epAI_actual = 4.5
epAI_pred = -4.5
epAI_ref = constants['ep_AI']

ref = mut.thermo.SimpleRepression(R=100, ep_r=ep_RA_ref,
                                 ka=constants['Ka'], ki=constants['Ki'],
                                 ep_ai=epAI_ref, n_sites=constants['n_sites'], 
                                 effector_conc=50)

ref_smooth = mut.thermo.SimpleRepression(R=100, ep_r=ep_RA_ref,
                                 ka=constants['Ka'], ki=constants['Ki'],
                                 ep_ai=epAI_ref, n_sites=constants['n_sites'], 
                                 effector_conc=c_range)

data = mut.thermo.SimpleRepression(R=100, ep_r=epRA_actual,
                                  ka=constants['Ka'], ki=constants['Ki'],
                                  ep_ai=epAI_actual, n_sites=constants['n_sites'],
                                  effector_conc=c_range_data)

prediction = mut.thermo.SimpleRepression(R=100, ep_r=ep_RA_pred,
                                  ka=constants['Ka'], ki=constants['Ki'],
                                  ep_ai=epAI_pred, n_sites=constants['n_sites'],
                                  effector_conc=c_range)

prediction_data = mut.thermo.SimpleRepression(R=100, ep_r=ep_RA_pred,
                                  ka=constants['Ka'], ki=constants['Ki'],
                                  ep_ai=epAI_pred, n_sites=constants['n_sites'],
                                  effector_conc=c_range_data)

true_deltaF = -ref.bohr_parameter() + data.bohr_parameter() 
predicted_deltaF = -ref.bohr_parameter() + prediction_data.bohr_parameter()

# Define the range of delta 
fig, ax = plt.subplots(2, 2, figsize=(5, 4.5))

# Plot the fake data
_ = ax[0, 0].plot(c_range / constants['Ka'], ref_smooth.fold_change(), color=colors[0], label='reference')
_ = ax[0, 0].plot(c_range / constants['Ka'], prediction.fold_change(), '-', color='firebrick', label='predicted')
_ = ax[0, 0].plot(c_range_data / constants['Ka'], data.fold_change(), '.', color='dodgerblue', label='observed')

# Label the first plot
_ = ax[0, 0].set_xscale('log')
_ = ax[0, 0].set_ylim([0, 1.1])
_ = ax[0, 0].set_xlabel('$c / K_A$')
_ = ax[0, 0].set_ylabel('fold-change')
_ = ax[0, 0].set_title('theory-experiment failure', backgroundcolor='#f1f2f6', y=1.08)
_leg = ax[0, 0].legend(fontsize=8, handlelength=1, bbox_to_anchor=(0.6, 0.5))

# Plot the coupling modes
_ = ax[0, 1].plot(deltaF_meas, deltaF_linear, color='firebrick')
_ = ax[1, 0].plot(deltaF_meas, deltaF_mult, color='firebrick')
_ = ax[1, 1].plot(deltaF_meas, deltaF_comp, color='firebrick')

titles = [ 'additive coupling\n$\Delta F_{true} = \Delta F_{pred} + \epsilon$', 
          'multiplicative coupling\n$\Delta F_{true} = \epsilon\Delta F_{pred}$', 
          'complicated coupling\n$\Delta F_{true} = \Delta F_{pred}(\epsilon)$']
for i, a in enumerate(ax.ravel()[1:]):
    # Plot agreement curve
    _ = a.plot(deltaF_range, deltaF_range, lw=1, ls=':')
    _ = a.set_xlim([-8, 8])
    _ = a.set_ylim([-8, 8])
    
    # Apply labels. 
    _ = a.set_xlabel('$\Delta F_{pred}$ [$k_BT$]')
    _ = a.set_ylabel('$\Delta F_{true}$ [$k_BT$]')
    
    # Add titles
    _ = a.set_title(titles[i], fontsize=8, backgroundcolor='#f1f2f6', y=1.08)

# Plot the observed data over the coupling modes
for a in ax.ravel()[1:]:
    _ = a.plot(predicted_deltaF, true_deltaF, '.', color='dodgerblue', label='observed', zorder=100)
    
# Add panel labels.
labels = ['(a)', '(b)', '(c)', '(d)']
for i, a in enumerate(ax.ravel()):
    a.text(-0.3, 1.2, labels[i], transform=a.transAxes, fontsize=8)
    
plt.tight_layout()
plt.savefig('./coupling_modes.pdf')



plt.figure()
plt.plot(predicted_deltaF, predicted_deltaF - true_deltaF)