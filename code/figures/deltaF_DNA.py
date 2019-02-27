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
DNA = data[(data['class']=='DNA') & (data['operator']=='O2')].copy()


# Find the inference statsitics for the  DNA binding energy
stats = pd.read_csv('../../data/csv/DNA_binding_energy_summary.csv')
stats = stats[(stats['operator']=='O2') & 
              (stats['repressors']==260) & (stats['parameter']=='ep_RA')][
              ['mutant', 'median', 'hpd_min','hpd_max']].copy()
ep_r = [constants[op] for op in DNA['operator']]
wt_bohr = -mut.thermo.SimpleRepression(DNA['repressors'], ep_r, ka=constants['Ka'],
                                    ki=constants['Ki'], ep_ai=constants['ep_AI'],
                                    effector_conc=DNA['IPTGuM']).bohr_parameter()
                                    
def delta_F_lims(bohr_min, bohr_max, bohr_ref):
    diff_min = [np.min([np.diff([ref, _min])[0],
        np.diff([ref, _max])[0]]) for ref, _min, _max  in zip(bohr_ref, bohr_min, bohr_max)]
    diff_max = [np.max([np.diff([ref, _min])[0],
                np.diff([ref, _max])[0]]) for ref, _min, _max  in zip(bohr_ref, 
                                                         bohr_min, bohr_max)]
    return [diff_min, diff_max]

_min, _max = delta_F_lims(DNA['bohr_min'].values, 
            DNA['bohr_max'].values, wt_bohr)
DNA['delta_F'] = wt_bohr - DNA['bohr_median']
DNA['delta_F_min'] = _min
DNA['delta_F_max'] = _max

# Define the marker shapes for plotting
rep_glyphs = {60:'s', 124:'d', 
              260:'o', 1220:'^'}


# Instantiate the figure canvas. 
# ############################
# PLOTTING
# ############################
fig, ax = plt.subplots(4, 3, figsize=(7, 5), sharex=True, sharey=True)
mut_idx = {'Q21M': 0, 'Q21A':1, 'Y20I':2}
rep_idx = {60:0, 124:1, 260:2, 1220:3}

# ############################
# DATA AND THEORY
# ###########################
for g, d in DNA.groupby(['mutant', 'repressors']):
    _ax = ax[rep_idx[g[1]], mut_idx[g[0]]]
    if g[1] == 260:
        face = 'w' 
    else:
        face = colors[g[0]] 

    _ = _ax.plot(d['IPTGuM'], d['delta_F'], marker=rep_glyphs[g[1]], 
                color = colors[g[0]], linestyle='none', ms=4, 
                markerfacecolor=face)
    _ = _ax.plot(d['IPTGuM'], d['delta_bohr_median'], marker=rep_glyphs[g[1]], 
    
                color = colors[g[0]], linestyle='none', ms=4, 
                markerfacecolor=face)

    _ = _ax.vlines(d['IPTGuM'], -d['delta_F_max'], -d['delta_F_min'], 
                   lw=0.75, color=colors[g[0]])

    median, hpd_min, hpd_max = -(constants['O2'] - 
    stats[stats['mutant']==g[0]][['median', 'hpd_min', 'hpd_max']].values[0])
    _ = _ax.hlines(median, -1, 1E4, color=colors[g[0]], lw=0.75)
    _ = _ax.fill_between([-1, 1E4], hpd_min, hpd_max, 
                        color=colors[g[0]], alpha=0.5)

 ################################
# FORMATTING
# ###############################
for a in ax.ravel():
    a.set_xscale('symlog')
    a.xaxis.set_tick_params(labelsize=9)
    a.yaxis.set_tick_params(labelsize=9)
    a.set_xlim([-1, 1E4])
    a.set_ylim([-8, 8])
    a.hlines(0, -1, 1E4, lw=0.75, color='slategray', linestyle='--',
    zorder=1)

# ##############################
# LABELING
# ###############################
for i in range(3):
    ax[3, i].set_xlabel('IPTG [µM]')
for i in range(4):
    ax[i, 0].set_ylabel('$\Delta F$ [$k_BT$]')
for m, ind in mut_idx.items():
    ax[0, ind].set_title(m, backgroundcolor=pboc['pale_yellow'], y=1.03,
    fontsize=9)
for r, ind in rep_idx.items():
    ax[ind, 0].text(-0.5, 0.52, str(int(r)), fontsize=9, rotation='vertical',
                    backgroundcolor=pboc['pale_yellow'], 
                    transform=ax[ind,0].transAxes)

plt.subplots_adjust(wspace=0.03, hspace=0.07)
fig.text(0.4, 0.98, 'DNA binding mutation', fontsize=10, backgroundcolor='#E3DBD0')
fig.text(-0.07, 0.6, 'repressors per cell', fontsize=10, backgroundcolor='#E3DBD0',
rotation='vertical')
# plt.tight_layout()
plt.savefig('../../figures/Fig3_DNA_deltaF.pdf', bbox_inches='tight')
