# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.thermo
import mut.stats
import mut.viz
colors = mut.viz.color_selector('mut')
pboc = mut.viz.color_selector('pboc')
mut.viz.plotting_style()
constants = mut.thermo.load_constants()

# Load the experimental data
data = pd.read_csv('../../data/csv/summarized_data.csv')
data = data[(data['class']=='IND')].copy()

# Load the inference data. 
kaki_samples = pd.read_csv('../../data/csv/Fig3_O2_KaKi_only_samples.csv')
kaki_stats = pd.read_csv('../../data/csv/Fig3_O2_KaKi_only_stats.csv')
kaki_allo_samples = pd.read_csv('../../data/csv/Fig3_O2_KaKi_epAI_samples.csv')
kaki_allo_stats = pd.read_csv('../../data/csv/Fig3_O2_KaKi_epAI_stats.csv')

# Define the plotting constants
c_range = np.logspace(-2, 4, 200)

# Instantiate the figure. 
fig, ax = plt.subplots(2, 4, figsize=(7,  4))

# Format the axes 
ax[0, 0].axis('off')
ax[1, 0].axis('off')
for a in ax.ravel():
    a.set_ylim([-0.05, 1.2])
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)

# Correct scaling. 
ax[0, 1].set_xscale('log')
ax[1, 1].set_xscale('log')

# Add appropriate labels
for i in range(2):
    ax[i, 1].set_xlabel('IPTG [M]', fontsize=8)
    ax[i, 1].set_ylabel('fold-change', fontsize=8)
    ax[i, 2].set_xlabel('Bohr parameter [$k_BT$]', fontsize=8)
    ax[i, 2].set_ylabel('fold-change', fontsize=8)
    ax[i, 1].set_xlim([1E-8, 1E-2])
    ax[i, 2].set_xlim([-6, 6])
    
# Compute and plot the collapse curves.
bohr_range = np.linspace(-6, 6, 200)
collapse = (1 + np.exp(-bohr_range))**-1
ax[0, 2].plot(bohr_range, collapse, lw=1, color='k')
ax[1, 2].plot(bohr_range, collapse, lw=1, color='k')

# Group and plt the data
grouped = data[data['operator']=='O2'].groupby(['mutant'])
for g, d in grouped:
    
    # Extract the statistics
    Ka = kaki_stats[kaki_stats['parameter']=='Ka.{}'.format(g)]['mode'].values[0]
    Ki = kaki_stats[kaki_stats['parameter']=='Ki.{}'.format(g)]['mode'].values[0]
    Ka_samples = kaki_samples['Ka.{}'.format(g)].values
    Ki_samples = kaki_samples['Ki.{}'.format(g)].values
    
    if (g == 'Q294K') | (g == 'Q294R'):
        _a = ax[1, 1]  
        _ac = ax[1, 2]
        plotcred=False
        ls= '--'
        
        # Extract the allosteric fit statistics.  
        Ka_allo = kaki_allo_stats[kaki_allo_stats['parameter']=='Ka.{}'.format(g)]['mode'].values[0]
        Ki_allo = kaki_allo_stats[kaki_allo_stats['parameter']=='Ki.{}'.format(g)]['mode'].values[0]
        epAI = kaki_allo_stats[kaki_allo_stats['parameter']=='ep_AI.{}'.format(g)]['mode'].values[0]
        Ka_allo_samples = kaki_allo_samples['Ka.{}'.format(g)].values
        Ki_allo_samples = kaki_allo_samples['Ki.{}'.format(g)].values
        epAI_samples = kaki_allo_samples['ep_AI.{}'.format(g)].values
       
        # Compute the modified bohr parameter. 
        bohr_allo = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=constants['O2'],
                                               ka=Ka_allo, ki=Ki_allo, ep_ai=epAI, 
                                               effector_conc=d['IPTGuM'], 
                                               n_sites=constants['n_sites'],
                                               n_ns=constants['Nns'],
                                               ).bohr_parameter()
        _ac.errorbar(bohr_allo, d['mean'], d['sem'], fmt='o', markerfacecolor='w', 
                          markeredgecolor=colors[g], markeredgewidth=1, markersize=3, 
                          linestyle='none', lw=1, color=colors[g], zorder=1000,
                        label='__nolegend__')
        bohr_cred = np.zeros([2, len(d['IPTGuM'])])
        for i, c in enumerate(d['IPTGuM']):
            _bohr = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=constants['O2'],
                                                   ka=Ka_allo_samples, ki=Ki_allo_samples, 
                                                    ep_ai=epAI_samples, effector_conc=c,
                                                   n_sites=constants['n_sites'],
                                                   n_ns=constants['Nns']).bohr_parameter()
            bohr_cred[:,  i] = mut.stats.compute_hpd(_bohr, mass_frac=0.95)
        _ac.hlines(d['mean'], bohr_cred[0, :], bohr_cred[1, :], lw=1, color=colors[g],
                  label='__nolegend__')
        
        # Compute the best-fit lines
        fit = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=constants['O2'],
                                     ka=Ka_allo, ki=Ki_allo, ep_ai=epAI,
                                     effector_conc=c_range, n_sites=constants['n_sites'], 
                                     n_ns=constants['Nns']).fold_change()    
        
        fc_cred = np.zeros([2, len(c_range)])
        for i, c in enumerate(c_range):
            # Compute the credible region. 
            _fit = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=constants['O2'],
                                     ka=Ka_allo_samples, ki=Ki_allo_samples, ep_ai=epAI_samples,
                                     effector_conc=c, n_sites=constants['n_sites'], 
                                     n_ns=constants['Nns']).fold_change()    
            fc_cred[:, i] = mut.stats.compute_hpd(_fit, mass_frac=0.95)
            
        _a.plot(c_range / 1E6, fit, color=colors[g])
        _a.fill_between(c_range / 1E6, fc_cred[0, :], fc_cred[1, :], alpha=0.2, color=colors[g])

    else: 
        ls = '-'
        _a = ax[0, 1]
        _ac = ax[0, 2]
        plotcred = True
    
    # Compute the fits. 
    fit = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=constants['O2'],
                                     ka=Ka, ki=Ki, ep_ai=constants['ep_AI'],
                                     effector_conc=c_range, n_sites=constants['n_sites'],
                                     n_ns=constants['Nns']).fold_change()
    
    
    # Plot the titration data. 
    _ = _a.errorbar(d['IPTGuM'] / 1E6, d['mean'], d['sem'], fmt='o', ms=3, linestyle='none',
                         linewidth=1, label=g.upper(), color=colors[g.upper()])

    # Plot the fits. 
    _ =_a.plot(c_range/1E6, fit, ls=ls, lw=1, color=colors[g])
    
    # Compute the credible regions
    if plotcred:
        fc_cred = np.zeros([2, len(c_range)])
        for i, c in enumerate(c_range):
            # Compute the credible region. 
            _fit = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=constants['O2'],
                                     ka=Ka_samples, ki=Ki_samples, ep_ai=constants['ep_AI'],
                                     effector_conc=c, n_sites=constants['n_sites'], 
                                     n_ns=constants['Nns']).fold_change()    
            fc_cred[:, i] = mut.stats.compute_hpd(_fit, mass_frac=0.95)
        _a.fill_between(c_range / 1E6, fc_cred[0, :], fc_cred[1, :], color=colors[g], alpha=0.2)
    # Compute the bohr parameter 
    bohr = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=constants['O2'],
                                ka=Ka, ki=Ki, ep_ai=constants['ep_AI'], 
                                effector_conc=d['IPTGuM'], n_sites=constants['n_sites'],
                                n_ns=constants['Nns']).bohr_parameter()   
    bohr_cred = np.zeros([2, len(d['IPTGuM'])])
    for i, c in enumerate(d['IPTGuM']):
            _bohr = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=constants['O2'],
                                                   ka=Ka_samples, ki=Ki_samples, 
                                                    ep_ai=constants['ep_AI'], effector_conc=c,
                                                   n_sites=constants['n_sites'],
                                                   n_ns=constants['Nns']).bohr_parameter()
            bohr_cred[:,  i] = mut.stats.compute_hpd(_bohr, mass_frac=0.95)
    _ac.errorbar(bohr, d['mean'], d['sem'], color=colors[g], linestyle='none', lw=1, fmt='o', 
                markersize=3, label='__nolegend__')
    
    _ac.hlines(d['mean'], bohr_cred[0, :], bohr_cred[1, :], lw=1, color=colors[g],
              label='__nolegend__')
    
# Format and add the legends. 
ax[1, 1].plot([], [], 'k--', label='no biophysical\nepistasis', lw=1)
ax[1, 1].plot([], [], 'k-', label='biophysical\nepistasis', lw=1)
ax[1, 2].plot([], [], 'o', markerfacecolor='w', markeredgecolor='k', markeredgewidth=1, markersize=3, 
              label='biophysical\nepistasis', lw=1)
ax[1, 2].plot([], [], 'o', markerfacecolor='k', markeredgecolor='k', markeredgewidth=1, markersize=3,
              label='no biophysical\nepistasis', lw=1)
ax[0, 1].legend(fontsize=7, handletextpad=0.1)
ax[1, 1].legend(bbox_to_anchor=(0.6, 0.35), fontsize=7, handletextpad=0.1)
ax[1, 2].legend(loc='upper left', fontsize=7, handletextpad=0.1)
plt.tight_layout()
plt.savefig('Fig3_IND_muts.svg', bbox_inches='tight')
