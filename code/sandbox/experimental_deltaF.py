# -*- coding: utf-8 -*-
import sys
import seaborn as sns
sys.path.insert(0, '../../')
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import mut.viz
import imp
imp.reload(mut.viz)
colors = mut.viz.personal_style()
constants = mut.thermo.load_constants()
constants['O3'] = -9.3

# Load the MWC induction data. 
data = pd.read_csv('../../data/csv/RazoMejia2018_data.csv')
data = data[data['repressors'] > 0]


# Define the reference architechture.
c_0 = 50 # in µM
ref_arch = mut.thermo.SimpleRepression(R=constants['RBS1'], ep_r=constants['O2'], ka=constants['Ka'],
                                      ki=constants['Ki'], ep_ai=constants['ep_AI'],
                                      effector_conc=c_0, n_sites=constants['n_sites'])

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
 
# Compute the theoretical curves. 
c_range = np.logspace(-2, 4, 200)
R_range = np.logspace(0, 4, 200)
epRA_range = np.linspace(-20, -5, 200)
ref_params = dict(effector_conc=c_0, Ka=constants['Ka'], Ki=constants['Ki'],
                  ep_AI=constants['ep_AI'], n_sites=constants['n_sites'],
                  ep_RA=constants['O2'], R=constants['RBS1'])
c_dep = dbohr(ref_params, dict(effector_conc=c_range, Ka=constants['Ka'], Ki=constants['Ki'],
                              ep_AI=constants['ep_AI'], n_sites=constants['n_sites'],
                              ep_RA =-13.9, R=constants['RBS1']))
ep_dep = dbohr(ref_params, dict(effector_conc=c_0, Ka=constants['Ka'], Ki=constants['Ki'],
                              ep_AI=constants['ep_AI'], n_sites=constants['n_sites'],
                              ep_RA =epRA_range, R=constants['RBS1']))
R_dep = dbohr(ref_params, dict(effector_conc=c_0, Ka=constants['Ka'], Ki=constants['Ki'],
                              ep_AI=constants['ep_AI'], n_sites=constants['n_sites'],
                              ep_RA =-13.9, R=R_range))


# Iterate through the data and compute the empirical deltaF
arch = mut.thermo.SimpleRepression(R=data['repressors']*2, ep_r=data['binding_energy'],
                                  ka=constants['Ka'], ki=constants['Ki'], ep_ai=constants['ep_AI'],
                                  n_sites=constants['n_sites'], n_ns=constants['Nns'],
                                  effector_conc=data['IPTG_uM'])

data['bohr_parameter'] = arch.bohr_parameter()
data['computed_F'] = np.log((1 / data['fold_change_A']) - 1)
data['empirical_delta_F'] = -ref_arch.bohr_parameter() - data['computed_F']

# Separate data sets for plotting and compute summary stats
c_data = data[(data['operator']=='O2') & (data['rbs']=='RBS1')]
c_data = c_data.groupby(['IPTG_uM'])['empirical_delta_F'].agg(('mean', 'sem')).reset_index()
ep_data = data[(data['rbs']=='RBS1') & (data['IPTG_uM']==c_0)]
ep_data = ep_data.groupby(['binding_energy'])['empirical_delta_F'].agg(('mean', 'sem')).reset_index()
r_data = data[(data['operator']=='O2') & (data['IPTG_uM']==c_0)]
r_data = r_data.groupby(['repressors'])['empirical_delta_F'].agg(('mean', 'sem')).reset_index()

# Set up the figure canvas. 
fig, ax = plt.subplots(1, 3, figsize=(7, 2.75))


# Plot the theoretical curves. 
_ = ax[0].plot(c_range, c_dep)
_ = ax[1].plot(epRA_range, ep_dep)
_ = ax[2].plot(R_range, R_dep)

# Plot a subset of the data. 
_ = ax[0].errorbar(c_data['IPTG_uM'], c_data['mean'], c_data['sem'], marker='.',
                  capsize=1, linestyle='none', color='firebrick', lw=1)
_ = ax[1].errorbar(ep_data['binding_energy'], ep_data['mean'], ep_data['sem'], marker='.',
                  capsize=1, linestyle='none', color='firebrick', lw=1)
_ = ax[2].errorbar(r_data['repressors'] * 2, r_data['mean'], r_data['sem'], marker='.',
                  capsize=1, linestyle='none', color='firebrick', lw=1)

# Add appropriate labeling
ax[0].set_xlabel('IPTG [µM]')
ax[0].set_ylabel(r'$\Delta F$ [$k_BT$]')
ax[1].set_ylabel(r'$\Delta F$ [$k_BT$]')
ax[2].set_ylabel(r'$\Delta F$ [$k_BT$]')
ax[1].set_xlabel('DNA binding\nenergy [$k_BT$]')
ax[2].set_xlabel('repressors per cell')

# Add panel labels. 
ax[0].text(-0.2, 1.15, '(a)', transform=ax[0].transAxes,
          fontsize=8)
ax[1].text(-0.2, 1.15, '(b)', transform=ax[1].transAxes,
          fontsize=8)
ax[2].text(-0.2, 1.15, '(c)', transform=ax[2].transAxes,
          fontsize=8)

# Add appropriate titles. 
ax[0].set_title(r'$R=R_0\,;\, \Delta\varepsilon_{RA}=\Delta\varepsilon_{RA_0}$',
               backgroundcolor='#f1f2f6', y=1.08)
ax[1].set_title(r'$R=R_0\,;\,c=c_0$', backgroundcolor='#f1f2f6', y=1.08)
ax[2].set_title(r'$c=c_0\,;\, \Delta\varepsilon_{RA}=\Delta\varepsilon_{RA_0}$',
               backgroundcolor='#f1f2f3', y=1.08)

# Add appropriate axis formatting. 
ax[0].set_xscale('log')
ax[2].set_xscale('log')

plt.tight_layout()
plt.savefig('./data_deltaF.pdf', bbox_inches='tight')