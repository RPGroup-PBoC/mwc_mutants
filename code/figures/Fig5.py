# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.thermo
import mut.stats
import mut.viz
color = mut.viz.color_selector('mut')
pboc = mut.viz.color_selector('pboc')
mut.viz.plotting_style()
constants = mut.thermo.load_constants()

# Load the data. 
data = pd.read_csv('../../data/csv/summarized_data.csv')
data = data[data['class'] == 'DBL']

indiv_samps = pd.read_csv('../../data/csv/Fig2_O2_DBL_samples.csv')
indiv_stats = pd.read_csv('../../data/csv/Fig2_O2_DBL_stats.csv')
predicted = pd.read_csv('../../data/csv/Fig5_O2_DBL_predicted_deltaBohr.csv')


# Load the various statistics
epRA_stats = pd.read_csv('../../data/csv/Fig2_O2_DNA_binding_energy_stats.csv')
kaki_stats = pd.read_csv('../../data/csv/Fig3_O2_KaKi_only_stats.csv')
allo_stats = pd.read_csv('../../data/csv/Fig3_O2_KaKi_epAI_stats.csv')

# Instantiate the figure
fig, ax = plt.subplots(1, 2, figsize=(7,4))

# Assign the labels and format. 
for a in ax:
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)

ax[0].set_ylabel('fold-change', fontsize=8)
ax[0].set_xlabel('Bohr parameter [$k_BT$]', fontsize=8)
ax[1].set_ylabel(r'$\Delta F_{c\rightarrow\infty}$', fontsize=8)

# Compute and plot the master curve
bohr_range = np.linspace(-8, 8, 200)
theo = (1 + np.exp(-bohr_range))**-1
ax[0].plot(bohr_range, theo, 'k-', label='__nolegend__')

# Compute the wild-type bohr parameter for all quereied c
wt_bohr = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=constants['O2'], 
                                      ka=constants['Ka'], ki=constants['Ki'], 
                                      n_sites=constants['n_sites'], n_ns=constants['Nns'],
                                     ep_ai=constants['ep_AI'], effector_conc=data['IPTGuM'].unique()).bohr_parameter()

mut_loc = {m:i for i, m in enumerate(data['mutant'].unique())}
# For each mutant, compute the most-likely value of bohr parameter 
for g, d in data.groupby('mutant'):
    dna, ind = g.split('-')
    
    # Isolate the statistics for the prediction.
    epRA = epRA_stats[epRA_stats['parameter']=='ep_RA.{}.260'.format(dna)]['mode'].values[0]
    
    if ind == 'Q294K':
        ka = allo_stats[allo_stats['parameter']=='Ka.{}'.format(ind)]['mode'].values[0]
        ki = allo_stats[allo_stats['parameter']=='Ki.{}'.format(ind)]['mode'].values[0]
        epAI = allo_stats[allo_stats['parameter']=='ep_AI.{}'.format(ind)]['mode'].values[0]
    else:
        ka = kaki_stats[kaki_stats['parameter']=='Ka.{}'.format(ind)]['mode'].values[0]
        ki = kaki_stats[kaki_stats['parameter']=='Ki.{}'.format(ind)]['mode'].values[0]
        epAI = constants['ep_AI']
    
    # Compute the predicted sbohr parameter. 
    bohr = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=epRA, ka=ka, ki=ki, ep_ai=epAI,
                                      n_ns=constants['Nns'], n_sites=constants['n_sites'], 
                                      effector_conc=d['IPTGuM'].values).bohr_parameter()
    
    # Compute the difference from the wt for the prediction
    pred_dBohr = np.sum((bohr - wt_bohr)**2) / len(d)
    
    # Compute the measured bohr. 
    meas_ka = indiv_stats[indiv_stats['parameter']=='Ka.{}'.format(g)]['mode'].values[0]
    meas_ki = indiv_stats[indiv_stats['parameter']=='Ki.{}'.format(g)]['mode'].values[0]
    meas_epRA = indiv_stats[indiv_stats['parameter']=='ep_RA.{}'.format(g)]['mode'].values[0]
    meas_bohr = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=meas_epRA,
                                           ka=meas_ka, ki=meas_ki, ep_ai=constants['ep_AI'],
                                           n_ns=constants['Nns'], n_sites=constants['n_sites'],
                                           effector_conc=d['IPTGuM'].values).bohr_parameter()
    meas_dBohr = np.sum((meas_bohr - wt_bohr)**2) / len(d)  
    # Plot the delta bohr
    ax[1].plot(mut_loc[g] - 0.1, pred_dBohr, color=color[g], marker='o', ms=3)
    ax[1].plot(mut_loc[g] + 0.1, meas_dBohr, color=color[g], marker='D', ms=3)
    ax[0].errorbar(bohr, d['mean'], d['sem'], color=color[g], fmt='o', ms=3, lw=1, alpha=0.75, label=g)
ax[0].legend(fontsize=6) 
# ax[1].set_xticks(mut_loc.values())
# ax[1].set_xticklabels(mut_loc.keys())
plt.savefig('data_collapse.pdf', bbox_inches='tight')

