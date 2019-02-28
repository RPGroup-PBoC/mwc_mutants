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
ep_r = np.array([constants[op] for op in data['operator']])
wt_bohr = -mut.thermo.SimpleRepression(data['repressors'], ep_r, ka=constants['Ka'],
                                    ki=constants['Ki'], ep_ai=constants['ep_AI'],
                                    effector_conc=data['IPTGuM']).bohr_parameter()


def delta_F_lims(bohr_min, bohr_max, bohr_ref):
    diff_min = [np.min([np.diff([ref, _min])[0],
        np.diff([ref, _max])[0]]) for ref, _min, _max  in zip(bohr_ref, bohr_min, bohr_max)]
    diff_max = [np.max([np.diff([ref, _min])[0],
                np.diff([ref, _max])[0]]) for ref, _min, _max  in zip(bohr_ref, 
                                                         bohr_min, bohr_max)]
    return [diff_min, diff_max]

_min, _max = delta_F_lims(data['bohr_min'].values, 
            data['bohr_max'].values, wt_bohr)
data['bohr_ref'] = wt_bohr
data['delta_F'] = (wt_bohr -data['bohr_median']) 
data['delta_fc'] = wt_fc - data['fold_change_median']
data['delta_F_min'] = _min
data['delta_F_max'] = _max

# Find the inference statsitics for the  DNA binding energy
stats = pd.read_csv('../../data/csv/DNA_binding_energy_summary.csv')
stats = stats[(stats['operator']=='O2') & 
              (stats['repressors']==260) & (stats['parameter']=='ep_RA')][
              ['mutant', 'median', 'hpd_min','hpd_max']].copy()

DNA = data[(data['class']=='DNA') & (data['operator']=='O2')].copy()
wt = data[(data['class']=='WT') & (data['operator']=='O2')]

# Define the marker shapes for plotting
rep_glyphs = {60:'s', 124:'d', 
              260:'o', 1220:'^'}


# Instantiate the figure canvas. 
# ############################
# PLOTTING
# ############################
fig, ax = plt.subplots(3, 4, figsize=(6.5, 3), sharex=True, sharey=True)
mut_idx = {'Q21M': 0, 'Q21A':1, 'Y20I':2}
rep_idx = {60:3, 124:2, 260:1, 1220:0}

# ############################
# DATA AND THEORY
# ###########################
for g, d in DNA.groupby(['mutant', 'repressors']):
    _ax = ax[mut_idx[g[0]], rep_idx[g[1]]]
    if g[1] == 260:
        face = 'w' 
    else:
        face = colors[g[0]] 

    _ = _ax.plot(d['IPTGuM'], d['delta_F'], marker='o', 
                color = colors[g[0]], linestyle='none', ms=4, 
                markerfacecolor=face)
    _ = _ax.plot(wt['IPTGuM'], wt['delta_F'], 'k,')
     
    _ = _ax.vlines(d['IPTGuM'], -d['delta_F_max'], -d['delta_F_min'], 
                   lw=0.75, color=colors[g[0]])


    median, hpd_min, hpd_max = -(constants['O2'] - 
    stats[stats['mutant']==g[0]][['median', 'hpd_min', 'hpd_max']].values[0])
    _ = _ax.hlines(median, -0.1, 1E4, color=colors[g[0]], lw=0.75)
    _ = _ax.fill_between([-0.1, 1E4], hpd_min, hpd_max, 
                        color=colors[g[0]], alpha=0.5)

 ################################
# FORMATTING
# ###############################
for a in ax.ravel():
    a.set_xscale('symlog')
    a.xaxis.set_tick_params(labelsize=9)
    a.yaxis.set_tick_params(labelsize=9)
    a.set_xlim([-0.3, 1E4])
    a.set_ylim([-5, 8])
    a.hlines(0, -0.1, 1E4, lw=0.75, color='slategray', linestyle='--',
    zorder=1)

# ##############################
# LABELING
# ###############################
for i in range(4):
    ax[2, i].set_xlabel('IPTG [µM]')
    ax[2, i].set_xticks([0, 1E0, 1E2, 1E4])
for i in range(3):
    ax[i, 0].set_ylabel('$\Delta F$ [$k_BT$]')
for m, ind in mut_idx.items():
    ax[ind, 0].text(-0.7, 0.57, m, rotation='vertical', backgroundcolor=pboc['pale_yellow'],
    fontsize=9, transform=ax[ind, 0].transAxes)
for r, ind in rep_idx.items():
    ax[0, ind].set_title('$R = $' + str(int(r)), fontsize=9,y=1.03,
                    backgroundcolor=pboc['pale_yellow'])
                    

plt.subplots_adjust(wspace=0.04, hspace=0.07)
plt.savefig('../../figures/Fig3_DNA_deltaF.pdf', bbox_inches='tight')
