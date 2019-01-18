# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
constants = mut.thermo.load_constants()
colors = mut.viz.personal_style()

# Define the function for computing the delta borh
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
 

# Define the reference states
c_0 = 50 #in µM
R_0 = 260 
epRA_0 = -13.9

# Define smooth ranges for theoretical plotting
c_range = np.logspace(-4,5, 500)
ep_range = np.linspace(-20, -5, 500)
rep_range = np.logspace(0, 4, 500)
 
# ###################
# DATA 
# ####################
# Load the data and trim to experimental strains. 
data = pd.read_csv('../../data/csv/RazoMejia2018_data.csv')
data = data[(data['repressors'] > 0) & (data['IPTG_uM'] > 0)]

data['repressors'] *= 2

# Compute the measured bohr parameter
data['measured_F'] = np.log((1 / data['fold_change_A']) - 1)

# Compute the  EC50 for each strain
data['predicted_ec50'] = mut.thermo.SimpleRepression(R=data['repressors'], ep_r=data['binding_energy'],
                                                  effector_conc=data['IPTG_uM'], ka=constants['Ka'],
                                                  ki=constants['Ki'], ep_ai=constants['ep_AI']).ec50()

# Compute the three reference bohrs. 
data['conc_ref'] = -1 * mut.thermo.SimpleRepression(R=data['repressors'], ep_r=data['binding_energy'],
                                                  effector_conc=c_0, ka=constants['Ka'],
                                                  ki=constants['Ki'], ep_ai=constants['ep_AI']).bohr_parameter()
data['rep_ref'] = -1 * mut.thermo.SimpleRepression(R=R_0, ep_r=data['binding_energy'],
                                                  effector_conc=data['IPTG_uM'], ka=constants['Ka'],
                                                  ki=constants['Ki'], ep_ai=constants['ep_AI']).bohr_parameter()
data['ep_ref'] = -1 * mut.thermo.SimpleRepression(R=data['repressors'], ep_r=epRA_0,
                                                  effector_conc=data['IPTG_uM'], ka=constants['Ka'],
                                                  ki=constants['Ki'], ep_ai=constants['ep_AI']).bohr_parameter()

# Compute the three empirical delta Fs
data['emp_conc_deltaF'] = data['conc_ref'] - data['measured_F']
data['emp_rep_deltaF'] = data['rep_ref'] - data['measured_F']
data['emp_ep_deltaF'] = data['ep_ref'] - data['measured_F']

# Compute the summaries. 
grouped_c = data.groupby('IPTG_uM')['emp_conc_deltaF'].agg(('mean', 'median', 'std')).reset_index()
grouped_ep = data.groupby('binding_energy')['emp_ep_deltaF'].agg(('mean', 'median', 'std')).reset_index()
grouped_rep = data.groupby('repressors')['emp_rep_deltaF'].agg(('mean', 'median', 'std')).reset_index()

# ####################
# THEORY
# ####################
# Compute the theoretical curves 
ref_params = {'R':260, 'ep_RA':-13.9, 'effector_conc':c_0,
                 'Ka':constants['Ka'], 'Ki':constants['Ki'],
                 'ep_AI': constants['ep_AI'], 'n_sites':constants['n_sites']}
per_params = {'R':260, 'ep_RA':-13.9, 'effector_conc':c_range,
                 'Ka':constants['Ka'], 'Ki':constants['Ki'],
                 'ep_AI': constants['ep_AI'], 'n_sites':constants['n_sites']}
delta_bohr_conc = dbohr(ref_params, per_params)
ref_params.update({'R':R_0, 'effector_conc': 50})
per_params.update({'R': rep_range, 'effector_conc':50})
delta_bohr_rep = dbohr(ref_params, per_params)
ref_params.update({'R':260, 'effector_conc': 50, 'ep_RA':epRA_0})
per_params.update({'R': 260, 'effector_conc':50, 'ep_RA':ep_range })
delta_bohr_ep = dbohr(ref_params, per_params)

# ######################
# PLOTTING
# #######################
# Insantiate the figure. 
fig, ax = plt.subplots(1, 3, figsize=(6, 2.5), sharey=True)

# Plot the theoretical curves
_ = ax[0].plot(c_0 / c_range, delta_bohr_conc, '-', color='rebeccapurple', lw=1, label='theory')
_ = ax[1].plot(epRA_0 - ep_range, delta_bohr_ep, '-', color='rebeccapurple', lw=1, label='theory')
_ = ax[2].plot(R_0 / rep_range, delta_bohr_rep, '-', color='rebeccapurple', lw=1, label='theory')

# Plot the summarized data
_ = ax[0].errorbar(c_0 / grouped_c['IPTG_uM'], grouped_c['mean'], grouped_c['std'],
                    ms=5, lw=1, capsize=2,marker='.', linestyle='none', color=colors[3], label='mean ± std')
_ = ax[1].errorbar(epRA_0 - grouped_ep['binding_energy'], grouped_ep['mean'], grouped_ep['std'],
                    ms=5, lw=1, capsize=2,marker='.', linestyle='none', color=colors[3], label='mean ± std')
_ = ax[2].errorbar(R_0 / grouped_rep['repressors'], grouped_rep['mean'], grouped_rep['std'],
                    ms=5, lw=1, capsize=2,marker='.', linestyle='none', color=colors[3], label='mean ± std')


       
# #########$$#########################
# FORMATTING
# ####################################

# Apply formatting and labeling
ax[0].set_xscale('log')
ax[2].set_xscale('log')
ax[0].set_xlabel('$c_0 / c$')
ax[1].set_xlabel(r'$\Delta\varepsilon_{RA_0} - \Delta\varepsilon_{RA}$ [$k_BT$]')
ax[2].set_xlabel('$R_0 / R$')
ax[0].set_ylabel(r'$\Delta F$ [$k_BT$]')

# Add panel labeling.
fig.text(0.01, 0.95, '(a)', fontsize=8)
fig.text(0.37, 0.95, '(b)', fontsize=8)
fig.text(0.67, 0.95, '(c)', fontsize=8)

for a in ax[1:]:
    a.spines['left'].set_visible(False)
    a.yaxis.set_tick_params(color='w')

# Titles and ylabels
titles = [r'$R = R_0\,;\, \Delta\varepsilon_{RA}=\Delta\varepsilon_{RA_0}$',
          r'$R = R_0 \,;\,c=c_0$', r'$c=c_0\,;\,\Delta\varepsilon_{RA}=\Delta\varepsilon_{RA_0}$']
for i, a in enumerate(ax):
     a.set_title(titles[i], fontsize=8, backgroundcolor='#f1f2f6',
               y=1.08)
ax[0].legend(loc='upper right')
plt.subplots_adjust(wspace=0.25)
plt.tight_layout()
plt.savefig('./complete_data_deltaF.pdf', bbox_inches='tight')
    
    
    
    