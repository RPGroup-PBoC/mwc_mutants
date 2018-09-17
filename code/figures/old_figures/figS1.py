## -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from collections import OrderedDict
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
sys.path.insert(0, '../../')
import mut.viz
import mut.stats
import mut.thermo
mut.viz.plotting_style()
colors = mut.viz.color_selector('mut')
pboc_colors = mut.viz.color_selector('pboc')
colors['wt'] = 'k'

# Load the compiled data.
data = pd.read_csv('../../data/csv/compiled_data.csv', comment='#')
data = data[(data['class'] == 'DNA') | (data['class'] == 'WT')]
data = data[data['fold_change'] > 0]
df = data

 # we have DNA mutants, and we have inducer mutants
# single O2 operator, plot as a function of IPTG concentration
DNAmut_list = ['wt', 'Q21A', 'Q21M', 'Y20I']
mut_list = {'wt':0, 'Q21A':1, 'Q21M':2, 'Y20I':3, 'Q294V':4, 'F164T':5, 'Q294R':6}

df_mcmc = pd.read_csv('../../data/mcmc/NB_emcee_mutants_DNA_strict.csv')

stats = mut.stats.compute_statistics(df_mcmc, logprob_name='lnprobability')

# Inducer parameters
ea = -4.935
ei = 0.635

ka = np.exp(-ea)
ki = np.exp(-ei)

# Define the IPTG concentrations to evaluate
# IPTG = np.logspace(-8, -2, 75)* 1E6
IPTG = np.logspace(-7, -2, 100)* 1E6
IPTG_lin = np.array([0, 1E-7])


fig, [ax_head, ax] = plt.subplots(2, 1, figsize=(4.5, 6))
fig.subplots_adjust(hspace=0.05)
ax_head.set_aspect(0.3, anchor='S')
ax_head.axis('off')

df_group = df.groupby(['mutant', 'IPTGuM'])

##############################
# plot fold-change vs bohr
##############################
# data
for i, data in df_group:
    if i[0] in DNAmut_list:
        if i[0] == 'wt':
            ep_r = -13.9
        else:
            stat = stats[stats.parameter.str.contains(i[0])]
            ep_r = stat['mode']

        bohr = mut.thermo.SimpleRepression(R=260.0,
                            ep_r = ep_r,
                            effector_conc = i[1], ka = np.exp(-ea), ki = np.exp(-ei),
                            ep_ai=4.5, n_sites=2).bohr_parameter()

        ax.errorbar(bohr, data.fold_change.mean(),
                    yerr =data.fold_change.std()/np.sqrt(len(data.fold_change)),
                    color = colors[i[0]],
                     fmt='o', markersize='4')

# theory
for m in DNAmut_list:
    if m == 'wt':
        ep_r = -13.9
    else:
        stat = stats[stats.parameter.str.contains(m)]
        ep_r = stat['mode'].values[0]

    # plot theory - using fit binding energies to cover Bohr parameter
    # over entire range relevant to the data
    fold_change = mut.thermo.SimpleRepression(R=260.0, ep_r = ep_r,
                            effector_conc = IPTG, ka = ka, ki = ki,
                            ep_ai=4.5, n_sites=2).fold_change(pact=True)
    bohr = mut.thermo.SimpleRepression(R=260.0, ep_r = ep_r,
                            effector_conc = IPTG, ka = ka, ki = ki,
                            ep_ai=4.5, n_sites=2).bohr_parameter()

    ax.plot(bohr, fold_change, color = 'k', label = m)

    bohr_head = np.linspace(bohr.min(), bohr.max())
    ax_head.plot(bohr_head,1.8*np.ones(len(bohr_head))*mut_list[m],
                 color = colors[m])
    bohr_head_center = (bohr.min() + bohr.max())/2
    if m=='wt':
        ax_head.text(x=bohr_head_center-1.5,y=1.55*mut_list[m],s='wild-type',
                 color = 'k', bbox=dict(facecolor='white', linewidth=0))
    else:
        ax_head.text(x=bohr_head_center-0.5,y=1.55*mut_list[m],s=m,
                 color = 'k', bbox=dict(facecolor='white', linewidth=0))


ax.set_xlabel('Bohr parameter ($k_B T$)', fontsize=12)
ax.set_ylabel('fold-change', fontsize=12)
ax.set_xlim(-6.2,6.2)
ax_head.set_xlim(-6.2,6.2)
ax.set_ylim(0,1.1)

# # for plotting the legend (removing replicate labels)
# from collections import OrderedDict
# handles, labels = ax.get_legend_handles_labels()
# by_label = OrderedDict(zip(labels, handles))
# by_label['wild-type'] = by_label.pop('wt')
# plt.legend(by_label.values(), by_label.keys(), ncol=4,
#           loc='lower center')

# plt.tight_layout()
fig.savefig('output_figs/figS1.pdf', bbox_inches='tight')
