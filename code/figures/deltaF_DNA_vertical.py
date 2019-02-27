# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import mut.viz
import mut.thermo
mut.viz.plotting_style()
constants = mut.thermo.load_constants()
colors = mut.viz.color_selector('mut')
pboc = mut.viz.color_selector('pboc')

# Load the data and isolate to the DNA binding mutants. 
data = pd.read_csv('../../data/csv/empirical_F_statistics.csv')
DNA = data[data['class']=='DNA'].copy()

# Find the inference statsitics for the  DNA binding energy
stats = pd.read_csv('../../data/csv/DNA_binding_energy_summary.csv')
stats = stats[(stats['operator']=='O2') & 
              (stats['repressors']==260) & (stats['parameter']=='ep_RA')][
              ['mutant', 'median', 'hpd_min','hpd_max']].copy()

rep_glyphs = {60:'s', 124:'d', 
              260:'o', 1220:'^'}

# Compute the wild-type bohr parameter
wt_bohr = -mut.thermo.SimpleRepression(R=DNA['repressors'], ep_r=constants['O2'],
                                      ka=constants['Ka'], ki=constants['Ki'],
                                      ep_ai=constants['ep_AI'], 
                                  effector_conc=DNA['IPTGuM']).bohr_parameter()

# Compute the empirical F and find the maximum and minimum differences. 
DNA['wt_bohr'] = wt_bohr
DNA['delta_F'] = wt_bohr - DNA['bohr_median']
DNA['delta_F_max'] = wt_bohr - DNA['bohr_max'] 
DNA['delta_F_min'] = wt_bohr - DNA['bohr_min']

# Instantiate the figure canvas. 
# ############################
# PLOTTING
# ############################
fig, ax = plt.subplots(4, 1, figsize=(3.42, 5), sharex=True, sharey=True)
rep_idx = {60:0, 124:1, 260:2, 1220:3}

# ############################
# DATA AND THEORY
# ###########################
for g, d in DNA.groupby(['mutant', 'repressors']):
    _ax = ax[rep_idx[g[1]]]
    if g[1] == 260:
        face = 'w' 
    else:
        face = colors[g[0]] 

    _ = _ax.plot(d['IPTGuM'], d['delta_F'], marker=rep_glyphs[g[1]], 
                color = colors[g[0]], linestyle='none', ms=4, 
                markerfacecolor=face)
    _ = _ax.vlines(d['IPTGuM'], d['delta_F_min'], d['delta_F_max'], 
                   lw=0.75, color=colors[g[0]])

    median, hpd_min, hpd_max = -(constants['O2'] - 
    stats[stats['mutant']==g[0]][['median', 'hpd_min', 'hpd_max']].values[0])
    _ = _ax.hlines(median, 0, 1E4, color=colors[g[0]], lw=0.75)
    _ = _ax.fill_between([0, 1E4], hpd_min, hpd_max, 
                        color=colors[g[0]], alpha=0.5)


 ################################
# FORMATTING
# ###############################
for a in ax.ravel():
    a.set_xscale('symlog')
    a.xaxis.set_tick_params(labelsize=9)
    a.yaxis.set_tick_params(labelsize=9)
    a.set_xlim([-0.1, 1E4])
    a.set_ylim([-3, 8])
    a.hlines(0, -1, 1E4, lw=0.75, color='slategray', linestyle='--',
    zorder=1)

# ##############################
# LABELING
# ###############################
ax[3].set_xlabel('IPTG [µM]')
for i in range(4):
    ax[i].set_ylabel('$\Delta F$ [$k_BT$]')
for r, ind in rep_idx.items():
    ax[ind].text(-0.3, 0.52, str(int(r)), fontsize=9, rotation='vertical',
                    backgroundcolor=pboc['pale_yellow'], 
                    transform=ax[ind].transAxes)

plt.subplots_adjust(wspace=0.03, hspace=0.07)
fig.text(-0.18, 0.6, 'repressors per cell', fontsize=10, backgroundcolor='#E3DBD0',
rotation='vertical')
plt.savefig('../../figures/Fig3_DNA_deltaF_vertical.pdf', bbox_inches='tight')
