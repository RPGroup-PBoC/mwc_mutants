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


# Generate a data set. 
R_0 = 260
epRA_actual = -11
epRA_fit= -10.9
epRA_0 = -13.9
c_0 = 50 # in µM

c_range = np.logspace(-2, 4, 12)
c_range_smooth = np.logspace(-2, 4, 200)
data = mut.thermo.SimpleRepression(R=R_0, ep_r=epRA_actual,
                                  ka=constants['Ka'], ki=constants['Ki'],
                                  ep_ai=constants['ep_AI'], n_sites=constants['n_sites'],
                                  effector_conc=c_range).fold_change() 


# Compute the theoretical delta Bohr as a function of c
ref_params = {'R':R_0, 'ep_RA':epRA_0,
             'effector_conc':c_0,
             'Ka':constants['Ka'],
             'Ki':constants['Ki'],
             'ep_AI':constants['ep_AI'],
             'n_sites':constants['n_sites']}
fit_params = {'R':R_0, 'ep_RA':epRA_fit,
             'effector_conc':c_range,
              'ep_AI':constants['ep_AI'],
             'Ka':constants['Ka'], 'Ki':constants['Ki'],
             'n_sites':constants['n_sites']}
ref_bohr = mut.thermo.SimpleRepression(R=R_0, ep_r=epRA_0, effector_conc=c_0, ep_ai=constants['ep_AI'],
                           ka=constants['Ka'], ki=constants['Ki'], n_sites=constants['n_sites']).bohr_parameter()

emp_deltaF = -ref_bohr - np.log((1/data) - 1) 

delta_F = dbohr(ref_params, fit_params)

# Plot the data set

# compute the delta bohr based on 
fig, ax = plt.subplots(1, 3, figsize=(7, 3))
ax[0].plot(c_range, data, '.')
ax[1].plot(c_0 / c_range, delta_F, '-', color='rebeccapurple', alpha=0.5, lw=0.75)
ax[1].plot(c_0 / c_range, emp_deltaF, '.', color=colors[1])



ax[0].set_xscale('log')
ax[1].set_xscale('log')
ax[2].plot(delta_F, emp_deltaF, '.', color=colors[1])
ax[2].plot([-8, 8], [-8, 8], '-', color='rebeccapurple', alpha=0.5, lw=0.75)





# Real data. 




# # Load in the mutant data. 
data = pd.read_csv('../../data/csv/compiled_data.csv')
Q21M = data[(data['mutant'] == 'Q21A')]
Q21M_R260 = Q21M[(Q21M['repressors']==260) & (Q21M['operator']=='O2')]

# Compute the theoretical delta Bohr as a function of c
ref_params = {'R':R_0, 'ep_RA':epRA_0,
             'effector_conc':c_0,
             'Ka':constants['Ka'],
             'Ki':constants['Ki'],
             'ep_AI':constants['ep_AI'],
             'n_sites':constants['n_sites']}
fit_params = {'R':R_0, 'ep_RA':-11.1,
             'effector_conc':Q21M_R260['IPTGuM'],
              'ep_AI':constants['ep_AI'],
             'Ka':constants['Ka'], 'Ki':constants['Ki'],
             'n_sites':constants['n_sites']}

ref_bohr = mut.thermo.SimpleRepression(R=R_0, ep_r=epRA_0, effector_conc=c_0, ep_ai=constants['ep_AI'],
                           ka=constants['Ka'], ki=constants['Ki'], n_sites=constants['n_sites']).bohr_parameter()

emp_deltaF = -ref_bohr - np.log((1/Q21M_R260['fold_change']) - 1) 

delta_F = dbohr(ref_params, fit_params)


# compute the delta bohr based on 
fig, ax = plt.subplots(1, 3, figsize=(7, 3))
ax[0].plot(Q21M_R260['IPTGuM'], Q21M_R260['fold_change'], '.')
ax[1].plot(c_0 / Q21M_R260['IPTGuM'], delta_F, '-', color='rebeccapurple', alpha=0.5, lw=0.75)
ax[1].plot(c_0 / Q21M_R260['IPTGuM'], emp_deltaF, '.', color=colors[1])



ax[0].set_xscale('log')
ax[1].set_xscale('log')
ax[2].plot(delta_F, emp_deltaF, '.', color=colors[1])
ax[2].plot([-8, 8], [-8, 8], '-', color='rebeccapurple', alpha=0.5, lw=0.75)



