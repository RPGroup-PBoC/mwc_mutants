# For operating system interaction
import os
import glob
import datetime
import sys

# For scientific computing
import numpy as np
import pandas as pd
import scipy.special

# import mut class
import sys
sys.path.insert(0, '../../')
import mut.bayes
import mut.viz
mut.viz.plotting_style()
colors = mut.viz.color_selector('mut')
import mut.stats
import mut.thermo

# Useful plotting libraries
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns

from collections import OrderedDict


################################################################################
# Load all of the 2018 flow data.
################################################################################
flow_files = glob.glob('../processing/2018*flow*/output/*fold_change.csv')
dfs = [pd.read_csv(f, comment='#') for f in flow_files]
flow_data = pd.concat(dfs, axis=0)
flow_data = flow_data[(flow_data['fold_change'] >= -0.2) & (flow_data['fold_change'] <= 1.3)]
flow_data = flow_data[(flow_data.mutant != 'Q21M') | (flow_data.IPTGuM != 0.0)]

# Load the microscopy data
mic_files = glob.glob('../processing/2018*microscopy*/output/*fold_change.csv')
dfs = [pd.read_csv(f) for f in mic_files]
mic_data = pd.concat(dfs, axis=0)
mic_data['IPTGuM'] = 0.0

df = pd.concat([flow_data, mic_data], ignore_index=True)
df = df[df.mutant != 'Q294R']


# Now we remove the autofluorescence and delta values
df = df[(df.mutant != 'auto') & (df.mutant != 'delta') & (df.operator == 'O2')]

# Restart index
df = df.reset_index()

# only bother with R=260 data
df = df[df.strain == 'R260']

################################################################################
# load in the mcmc chains
################################################################################

df_mcmc = pd.read_csv('../../data/mcmc/NB_emcee_mutants_DNA_strict.csv')
stats = mut.stats.compute_statistics(df_mcmc, logprob_name='lnprobability')

################################################################################
# Plotting!
################################################################################

 # we have DNA mutants,
DNAmut_list = ['wt', 'Q21A', 'Q21M', 'Y20I']

# Inducer parameters
ea = -4.935
ei = 0.635

ka = np.exp(-ea)
ki = np.exp(-ei)

# Define the IPTG concentrations to evaluate
# IPTG = np.logspace(-8, -2, 75)* 1E6
IPTG = np.logspace(-7, -2, 100)* 1E6
IPTG_lin = np.array([0, 1E-7])

# Set the colors for the strains
# colors = sns.color_palette('colorblind', n_colors=7)
# colors[4] = sns.xkcd_palette(['dusty purple'])[0]
mut_list = {'wt':0, 'Q21A':1, 'Q21M':2, 'Y20I':3, 'Q294V':4, 'F164T':5, 'Q294R':6}


fig, [[leg,ax2_head], [ax1,ax2]] = plt.subplots(2, 2, figsize=(9, 6))
fig.subplots_adjust(hspace=0.05)
ax2_head.set_aspect(0.3, anchor='S')
ax2_head.axis('off')
# _.set_aspect(0.3, anchor='S')
# _.axis('off')

leg.set_aspect(0.3, anchor='S')
leg.axis('off')

df_group = df.groupby(['mutant', 'IPTGuM'])

##############################
# left plot - fold-change vs IPTG
##############################
# data
for i, data in df_group:
    if i[0] in DNAmut_list:
        ax1.errorbar(i[1]*10**-6, data.fold_change.mean(),
                     yerr =data.fold_change.std()/np.sqrt(len(data.fold_change)),
                     label=i[0], color = colors[i[0]],
                     fmt='o', markersize='4')
# theory
for m in DNAmut_list:
    stat = stats[stats.parameter.str.contains(m)]

    # plot average fold-change using fit binding energy
    fold_change = mut.thermo.SimpleRepression(R=260.0, ep_r=stat['mode'].values[0],
                            effector_conc = IPTG, ka = ka, ki = ki,
                            ep_ai=4.5, n_sites=2).fold_change(pact=True)
    ax1.plot(IPTG*10**-6, fold_change, color = colors[m])

        # Linear scale
    fold_change = mut.thermo.SimpleRepression(R=260.0, ep_r=stat['mode'].values[0],
                            effector_conc = IPTG_lin, ka = ka, ki = ki,
                            ep_ai=4.5, n_sites=2).fold_change(pact=True)
    ax1.plot(IPTG_lin, fold_change, color = colors[m],
            label=None, zorder=1, linestyle=':')

    # plot hpd bounds using fit binding energy
    fold_change_min = mut.thermo.SimpleRepression(R=260.0, ep_r=stat['hpd_min'].values[0],
                            effector_conc = IPTG, ka = np.exp(-ea), ki = np.exp(-ei),
                            ep_ai=4.5, n_sites=2).fold_change(pact=True)
    fold_change_max = mut.thermo.SimpleRepression(R=260.0, ep_r=stat['hpd_max'].values[0],
                            effector_conc = IPTG, ka = np.exp(-ea), ki = np.exp(-ei),
                            ep_ai=4.5, n_sites=2).fold_change(pact=True)

    ax1.fill_between(IPTG*10**-6, fold_change_min, fold_change_max,
                     alpha=0.3, color = colors[m])


# Set the sclae and labels.
ax1.set_xscale('symlog', linthreshx=1E-7, linscalex=0.5)
ax1.set_xlabel('IPTG (M)', fontsize=12)
ax1.set_ylabel('fold-change', fontsize=12)
ax1.set_ylim([-0.01, 1.1])
ax1.set_xlim([-5E-9, 1E-2])
ax1.tick_params(labelsize=12)



##############################
# right plot - fold-change vs bohr
##############################
# data
for i, data in df_group:
    if i[0] in DNAmut_list:
        stat = stats[stats.parameter.str.contains(i[0])]

        bohr = mut.thermo.SimpleRepression(R=260.0,
                            ep_r=stat['mode'],
                            effector_conc = i[1], ka = np.exp(-ea), ki = np.exp(-ei),
                            ep_ai=4.5, n_sites=2).bohr_parameter()

        ax2.errorbar(bohr, data.fold_change.mean(),
                    yerr =data.fold_change.std()/np.sqrt(len(data.fold_change)),
                    color = colors[i[0]],
                     fmt='o', markersize='4')




# theory
for m in DNAmut_list:
    stat = stats[stats.parameter.str.contains(m)]

    # plot theory - using fit binding energies to cover Bohr parameter
    # over entire range relevant to the data
    fold_change = mut.thermo.SimpleRepression(R=260.0, ep_r=stat['mode'].values[0],
                            effector_conc = IPTG, ka = ka, ki = ki,
                            ep_ai=4.5, n_sites=2).fold_change(pact=True)
    bohr = mut.thermo.SimpleRepression(R=260.0, ep_r=stat['mode'].values[0],
                            effector_conc = IPTG, ka = ka, ki = ki,
                            ep_ai=4.5, n_sites=2).bohr_parameter()

    ax2.plot(bohr, fold_change, color = 'k')

    bohr_head = np.linspace(bohr.min(), bohr.max())
    ax2_head.plot(bohr_head,1.8*np.ones(len(bohr_head))*mut_list[m],
                 color = colors[m])
    bohr_head_center = (bohr.min() + bohr.max())/2
    if m=='wt':
        ax2_head.text(x=bohr_head_center-1.5,y=1.55*mut_list[m],s='wild-type',
                 color = 'k', bbox=dict(facecolor='white', linewidth=0))
    else:
        ax2_head.text(x=bohr_head_center-0.5,y=1.55*mut_list[m],s=m,
                 color = 'k', bbox=dict(facecolor='white', linewidth=0))


ax2.set_xlabel('Bohr parameter ($k_B T$)', fontsize=12)
ax2.set_ylabel('fold-change', fontsize=12)
ax2.set_xlim(-6.2,6.2)
ax2_head.set_xlim(-6.2,6.2)
ax2.set_ylim(0,1.1)

# for plotting the legend (removing replicate labels)
handles, labels = ax1.get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
by_label['wild-type'] = by_label.pop('wt')
leg.legend(by_label.values(), by_label.keys(), ncol=4,
          loc='lower center')

# plt.tight_layout()
fig.savefig('output_figs/MWC_mutants_fig4.pdf', bbox_inches='tight')
