# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
colors = mut.viz.color_selector('mut')
pboc = mut.viz.color_selector('pboc')
rep_colors = {60:pboc['blue'], 124:pboc['purple'], 
              260:pboc['red'], 1220:pboc['green']}
constants = mut.thermo.load_constants()
mut.viz.plotting_style()

# ##############################################################################
# DATA LOADING / PRUNING
# ##############################################################################
data = pd.read_csv('../../data/csv/empirical_F_data_statistics.csv')
DNA = data[(data['class']=='WT') | (data['class']=='DNA')]
stats = pd.read_csv('../../data/csv/DNA_binding_energy_summary.csv')
stats = stats[(stats['repressors']==260) & (stats['operator']=='O2')]

# ##############################################################################
# FIGURE FORMATTING / INSTANTIATION
# ##############################################################################
fig, ax = plt.subplots(3, 4, figsize=(7, 3.5), sharey=True, sharex=True)
mut_axes = {'Y20I': 0, 'Q21A': 1, 'Q21M': 2}    
rep_axes = {60:0, 124: 1, 260:2, 1220:3}
for a in ax.ravel():
    a.set_ylim([-8, 8])
    a.set_xlim([-6, 4.5])
    a.set_xticks([-6, -4, -2, 0 , 2, 4])
    a.xaxis.set_tick_params(labelsize=7)
    a.yaxis.set_tick_params(labelsize=7)

for i in range(4):
    ax[-1, i].set_xlabel('$F^{(\mathrm{ref})}$ [$k_BT$]', fontsize=8)
for r, ind in rep_axes.items():
    ax[0, ind].set_title('$R = $'+  str(int(r)), 
                       backgroundcolor=pboc['pale_yellow'], y=1.08,fontsize=8)
for m, ind in mut_axes.items():
    ax[ind, 0].set_ylabel('$\Delta F$ [$k_BT$]', fontsize=8)
    ax[ind, 0].text(-0.45, 0.6, m, transform=ax[ind, 0].transAxes, 
                    rotation='vertical', backgroundcolor=pboc['pale_yellow'], 
                    fontsize=8)
    
# ##############################################################################
# EMPIRICAL ∆F PLOTS
# ############################################################################## 
for g, d in DNA[DNA['mutant']!='wt'].groupby(['mutant', 'repressors']):
    _ax = ax[mut_axes[g[0]], rep_axes[int(g[1])]] 
    if g[1] == 260:
        face = 'w'
    else:
        face=rep_colors[int(g[1])]
    dbohr = d[d['parameter']=='delta_bohr']
    wt_bohr = mut.thermo.SimpleRepression(R=g[1], ep_r=constants['O2'],
                                         ka=constants['Ka'], ki=constants['Ki'],
                                         n_sites=constants['n_sites'],
                                         ep_ai=constants['ep_AI'],
                                         effector_conc=dbohr['IPTGuM']).bohr_parameter()
    _ax.plot(wt_bohr, dbohr['median'], 'o', color=rep_colors[int(g[1])], ms=4,
            markerfacecolor=face, linestyle='none')
    _ax.vlines(wt_bohr, dbohr['hpd_min'], dbohr['hpd_max'], lw=1, color=rep_colors[int(g[1])])


# ################################################################
# PREDICTIONS 
# ################################################################
for g, d in stats[stats['mutant'] != 'wt'].groupby(['mutant']):
    _d =  d[d['parameter']=='ep_RA']
    for r, ind in rep_axes.items():
        _ax = ax[mut_axes[g], ind]
        bohr_range = np.linspace(_ax.get_xlim()[0],_ax.get_xlim()[1], 100)

        _ax.fill_between(bohr_range, constants['O2'] - _d['hpd_min'], constants['O2'] - _d['hpd_max'],
                        color=pboc['light_red'], zorder=1)
        
        _ax.set_xlim([bohr_range[0], bohr_range[-1]])
for a in ax.ravel():
    a.hlines(0, a.get_xlim()[0], a.get_xlim()[1], color='k', linestyle=':', lw=1)
plt.tight_layout()
plt.savefig('../../figures/Fig4_DNA_deltaF_dataref.pdf', bbox_inches='tight')