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

# Load the MWC induction data. 
data = pd.read_csv('../../data/csv/RazoMejia2018_data.csv')
data = data[data['repressors'] > 0]

# Define constants. 
c_range = np.logspace(-2, 5, 200)
R_range = np.logspace(0, 4, 200)
epRA_range = np.linspace(-18, -8, 200)

ref_arch = mut.thermo.SimpleRepression(R=constants['RBS1027'],
                                     ep_r=constants['O2'],
                                     ka=constants['Ka'], 
                                     ki=constants['Ki'],
                                     ep_ai=constants['ep_AI'],
                                     n_sites=constants['n_sites'],
                                     effector_conc=c_range)

# Set the function for computing delta bohr. 
def dbohr(ref_params, per_params):
    return np.log(ref_params['pact']/per_params['pact']) +\
           np.log(ref_params['R']/per_params['R']) -\
           (ref_params['ep_RA'] - per_params['ep_RA'])

# ###########################
# THEORETICAL CURVES
# ###########################
# Inducer concentration dependence
c0_pact = mut.thermo.MWC(effector_conc=ref_arch.ec50(),
                        ka=constants['Ka'], ki=constants['Ki'],
                        ep_ai=constants['ep_AI']).pact()
c0_bohr = mut.thermo.SimpleRepression(R=constants["RBS1027"], ep_r=constants['O2'],
                                     ka=constants['Ka'], ki=constants['Ki'],
                                     ep_ai=constants['ep_AI'], n_sites=constants['n_sites'],
                                     effector_conc=ref_arch.ec50()).bohr_parameter()
cstar_pact = mut.thermo.MWC(effector_conc=c_range,
                           ka=constants['Ka'], ki=constants['Ki'],
                           ep_ai=constants['ep_AI']).pact()
c0_theo = dbohr(dict(pact=c0_pact, R=constants['RBS1027'], ep_RA=constants['O2']),
                dict(pact=cstar_pact, R=constants['RBS1027'], ep_RA=constants['O2']))


# Repressor copy number dependence
R0_theo = dbohr(dict(pact=1, R=constants['RBS1027'], ep_RA=constants['O2']),
               dict(pact=1, R=R_range, ep_RA=constants['O2']))

# DNA binding energy dependence
epRA_theo = dbohr(dict(pact=1, R=constants['RBS1027'], ep_RA=constants['O2']),
                 dict(pact=1, R=constants['RBS1027'], ep_RA=epRA_range))


# #####################################
# HEATMAPS
# #####################################
# Induction profile
c0 = ref_arch.ec50()
r_mesh, c_mesh = np.meshgrid(R_range, cstar_pact)
c_r = dbohr(ref_params=dict(pact=c0_pact, R=260, ep_RA=-13.9),
            per_params=dict(pact=c_mesh, R=r_mesh, ep_RA=-13.9))
ep_mesh, c_mesh = np.meshgrid(epRA_range, cstar_pact)
c_ep = dbohr(ref_params=dict(pact=c0_pact, R=260, ep_RA=-13.9),
            per_params=dict(pact=c_mesh, R=260, ep_RA=ep_mesh))

#Repressor titration curve
r_pact = mut.thermo.MWC(effector_conc=c_range, ka=constants['Ka'],
                           ki=constants['Ki'], ep_ai=constants['ep_AI'],
                           n_sites=constants['n_sites']).pact()

r_mesh, ep_mesh = np.meshgrid(R_range, epRA_range)
r_c = dbohr(ref_params=dict(pact=r_pact, R=260, ep_RA=-13.9),
            per_params=dict(pact=c_mesh, R=r_mesh, ep_RA=-13.9))
ep_mesh, c_mesh = np.meshgrid(epRA_range, cstar_pact)
r_ep = dbohr(ref_params=dict(pact=1, R=260, ep_RA=-13.9),
            per_params=dict(pact=1, R=r_mesh, ep_RA=ep_mesh))

# #############################
# DATA 
# #############################
# Isolate the data sets
O2_R260_data = data[(data['operator']=='O2') & (data['rbs']=="RBS1027")]
O2_data = data[(data['operator']=='O2') & (data['IPTG_uM']==50)]
R260_data = data[(data['rbs']=='RBS1027') & (data['IPTG_uM']==50)]

# Compute the empirical ∆F
O2_R260_data['computed_F'] = -np.log((1/O2_R260_data['fold_change_A']) - 1)
O2_R260_data['deltaF'] = c0_bohr - O2_R260_data['computed_F']
R0_bohr = mut.thermo.SimpleRepression(R=constants["RBS1027"], ep_r=constants['O2'],
                                     ka=constants['Ka'], ki=constants['Ki'],
                                     ep_ai=constants['ep_AI'], n_sites=constants['n_sites'],
                                     effector_conc=O2_data['IPTG_uM']).bohr_parameter()
R260_bohr = mut.thermo.SimpleRepression(R=constants["RBS1027"], ep_r=constants['O2'],
                                     ka=constants['Ka'], ki=constants['Ki'],
                                     ep_ai=constants['ep_AI'], n_sites=constants['n_sites'],
                                     effector_conc=R260_data['IPTG_uM']).bohr_parameter()
O2_data['computed_F'] = -np.log((1/O2_data['fold_change_A']) - 1)
O2_data['deltaF'] = R0_bohr - O2_data['computed_F']
R260_data['computed_F'] = -np.log((1/R260_data['fold_change_A']) -1)
R260_data['deltaF'] = R260_bohr - R260_data['computed_F']


# ######################################
# CANVAS FORMATTING 
# ######################################
# Set up the figure canvas
fig, ax = plt.subplots(3, 3, figsize=(7, 6))
xlabels = ['IPTG [µM]', 'repressors per cell', 'DNA binding energy [$k_BT$]']
ax = ax.ravel()
for i, a in enumerate(ax.ravel()[:3]):
    if i != 2:
        a.set_xscale('log')
    a.set_ylabel(r'$\Delta F$ [$k_BT$]')
    a.set_xlabel(xlabels[i])

# Add appropriate labels. 
ax[0].set_title(r'$c^0 = [EC_{50}]$', backgroundcolor='#f1f2f6', y=1.08)
ax[1].set_title(r'$R^0 = 260$', backgroundcolor='#f1f2f6', y=1.08)
ax[2].set_title(r'$\Delta\varepsilon_{RA}^0 = -13.9\, k_BT$', backgroundcolor='#f1f2f6', y=1.08) 
_ = ax[0].plot(c_range, -c0_theo, label='theory')
_ = ax[1].plot(R_range, -R0_theo, label='theory')
_ = ax[2].plot(epRA_range, -epRA_theo, label='theory')

# #####################################
# DATA PLOTTING 
# #####################################
for g, d in O2_R260_data.groupby('IPTG_uM'):
    if g==0:
        label='data'
    else:
        label='__nolegend__'
    _ = ax[0].errorbar(g, d['deltaF'].mean(), d['deltaF'].std() / np.sqrt(len(d)), marker='.',
                      color=colors[2], lw=1, capsize=1, label=label, linestyle='none')

for g, d in O2_data.groupby(['repressors', 'IPTG_uM']):
    if g[0]==0:
        label='data'
    else:
        label='__nolegend__'
    _ = ax[1].errorbar(g[0] * 2, d['deltaF'].mean(), d['deltaF'].std() / np.sqrt(len(d)), marker='.',
                      color=colors[2], lw=1, capsize=1, label=label, linestyle='none')
for g, d in R260_data.groupby('operator'):
    if g=='O2':
        label='data'
    else:
        label='__nolegend__'
    _ = ax[2].errorbar(constants[g], d['deltaF'].mean(), d['deltaF'].std() / np.sqrt(len(d)), marker='.',
                      color=colors[2], lw=1, capsize=1, label=label, linestyle='none')
for a in ax.ravel()[:3]:
    a.legend()
    
# ####################################
# HEATMAP PLOTTING 
# ####################################
for a in ax[3:]:
    a.grid(False)
ax[3].imshow(c_r, cmap='Blues', interpolation='none', 
            aspect='auto')
ax[3].contour(c_r, levels=[0], linestyles=':', colors='w')
ax[3].set_xlabel('IPTG [µM]')
ax[3].set_ylabel('repressors / cell')

ax[4].imshow(r_c, cmap='Blues', interpolation='none',
            aspect='auto')
ax[4].contour(r_c, levels=[0], linestyles=':', colors='w')
ax[4].set_xlabel('repressors per cell')
ax[4].set_ylabel('IPTG[µM]')

ax[6].imshow(c_ep, cmap='Blues', interpolation='none', 
            aspect='auto')
ax[6].contour(c_ep, levels=[0], linestyles=':', colors='w')
ax[6].set_xlabel('IPTG [µM]')
ax[6].set_ylabel('DNA binding energy [$k_BT$]')

ax[7].imshow(r_ep, cmap='Blues', interpolation='none',
            aspect='auto')
ax[7].contour(r_ep, levels=[0], linestyles=':', colors='w')
ax[7].set_xlabel('repressors per cell')
ax[7].set_ylabel('DNA binding energy [$k_BT$]')
plt.tight_layout()
    
    

