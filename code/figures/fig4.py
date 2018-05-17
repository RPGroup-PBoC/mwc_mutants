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
FIG_NO = 4

# Load the compiled data.
data = pd.read_csv('../../data/csv/compiled_data.csv', comment='#')
data = data[(data['class'] == 'IND') | (data['class'] == 'WT')]
data = data[data['fold_change'] > 0]
df = data

# Load the fitting chains.
chains = pd.read_csv('../../../data/mcmc/NB_emcee_mutants_IND_strict.csv')
ind = np.argmax(chains['lnprobability'].values)
modes = pd.DataFrame(chains.drop(
    columns=['lnprobability', 'sigma', 'Unnamed: 0']).iloc[ind]).reset_index()
modes.columns = ['parameter', 'mode']

# determine parameter hpd
stats = mut.stats.compute_statistics(chains.drop(
    columns=['sigma', 'Unnamed: 0']), logprob_name='lnprobability')

# Define the architectural parameters.
Nns = 4.6E6
ka = 139E-6  # in M
ki = 0.53E-6  # in M
ep_r = {i.split('_')[0]: j for i, j in zip(modes['parameter'], modes['mode'])}
ep_r_hpd = [{i.split('_')[0]: j for i, j in zip(modes['parameter'], modes['mode'])},
            {i.split('_')[0]: j for i, j in zip(modes['parameter'], modes['mode'])}]
rep_list = data['repressors'].unique()
markers = ['o', 'D', 's', 'X']
marker_dict = dict(zip(rep_list,markers))

# Define the IPTG concentrations to evaluate
# IPTG = np.logspace(-8, -2, 75)* 1E6
IPTG = np.logspace(-7, -2, 100)* 1E6
IPTG_lin = np.array([0, 1E-7])* 1E6

INDmut_list = ['wt', 'Q294V', 'F164T']

# WT DNA parameters
ep_r = -13.9
# ea = -4.935
# ei = 0.635
# ka = np.exp(-ea)
# ki = np.exp(-ei)

# Define the IPTG concentrations to evaluate
IPTG = np.logspace(-7, -2, 100)* 1E6
IPTG_lin = np.array([0, 1E-7])

mut_list = {'wt':0, 'Q294V':1, 'F164T':2, 'Q294R':6}

fig = plt.figure(figsize=(5, 2.5))

ax = [[],[]]
ax[0] = plt.subplot(121)
ax[1] = plt.subplot(122)
df_group = df.groupby(['mutant', 'IPTGuM'])


##############################
# left plot - fold-change vs IPTG
##############################
# data
for i, data in df_group:
    if i[0] in INDmut_list:
        ax[0].errorbar(i[1]*10**-6, data.fold_change.mean(),
                     yerr =data.fold_change.std()/np.sqrt(len(data.fold_change)),
                     label=i[0], color = colors[i[0]],
                     fmt='o', markersize='4')
# theory
for m in INDmut_list:
    if m =='wt':
        ea = -4.935
        ei = 0.635
        stat_KA = np.exp(-ea)
        stat_KI = np.exp(-ei)
    else:
        stat = stats[stats.parameter.str.contains(m)]

        stat_KA = stat[stat.parameter.str.contains('Ka')]
        stat_KI = stat[stat.parameter.str.contains('Ki')]

    # plot average fold-change using fit binding energy
    fold_change = mut.thermo.SimpleRepression(R=260.0, ep_r=ep_r,
                            effector_conc = IPTG, ka = stat_KA['mode'].values[0],
                            ki = stat_KI['mode'].values[0],
                            ep_ai=4.5, n_sites=2).fold_change(pact=True)
    ax[0].plot(IPTG*1E-6, fold_change, color = colors[m])

    # Linear scale
    fold_change = mut.thermo.SimpleRepression(R=260.0, ep_r=ep_r,
                            effector_conc = IPTG_lin, ka = stat_KA['mode'].values[0],
                            ki = stat_KI['mode'].values[0],
                            ep_ai=4.5, n_sites=2).fold_change(pact=True)
    ax[0].plot(IPTG_lin, fold_change, color = colors[m],
            label=None, zorder=1, linestyle=':')

    # plot hpd bounds using fit binding energy
    fold_change_min = mut.thermo.SimpleRepression(R=260.0, ep_r=ep_r,
                            effector_conc = IPTG, ka = stat_KA['hpd_min'].values[0],
                            ki = stat_KI['hpd_min'].values[0],
                            ep_ai=4.5, n_sites=2).fold_change(pact=True)
    fold_change_max = mut.thermo.SimpleRepression(R=260.0,  ep_r=ep_r,
                            effector_conc = IPTG, ka = stat_KA['hpd_max'].values[0],
                            ki = stat_KI['hpd_max'].values[0],
                            ep_ai=4.5, n_sites=2).fold_change(pact=True)

    ax[0].fill_between(IPTG*10**-6, fold_change_min, fold_change_max,
                     alpha=0.3, color = colors[m])

# Set the sclae and labels.
ax[0].set_xscale('symlog', linthreshx=1E-7, linscalex=0.5)
ax[0].set_xlabel('IPTG (M)', fontsize=8)
ax[0].set_ylabel('fold-change', fontsize=8)
ax[0].set_ylim([-0.01, 1.1])
ax[0].set_xlim([-5E-9, 1E-2])
ax[0].tick_params(labelsize=8)


##############################
# right plot - contour K_A/K_I
##############################

#here's our data to plot, all normal Python lists
x = np.logspace(-5,-3)
y = np.logspace(-7,-5)

#setup the 2D grid with Numpy
x, y = np.meshgrid(x, y)
data = plt.cm.jet([x, y])
dynam = mut.thermo.SimpleRepression(R=260.0, ep_r=ep_r,
                        effector_conc = IPTG, ka = x,
                        ki = y,
                        ep_ai=4.5, n_sites=2).dynamic_range()

sc1 = ax[1].contourf(x, y, dynam, 10, vmin=0.0,vmax=1.0, cmap=cm.PuBu_r)

# data
for i, data in df_group:
    if i[0] in INDmut_list:
            if m =='wt':
                ea = -4.935
                ei = 0.635
                stat_KA = np.exp(-ea)
                stat_KI = np.exp(-ei)
            else:
                stat = stats[stats.parameter.str.contains(i[0])]

                stat_KA = stat[stat.parameter.str.contains('Ka')]
                stat_KI = stat[stat.parameter.str.contains('Ki')]

        dynam = mut.thermo.SimpleRepression(R=260.0,
                             ep_r=ep_r,
                            effector_conc = i[1], ka = stat_KA['mode'].values[0],
                            ki = stat_KI['mode'].values[0],
                            ep_ai=4.5, n_sites=2).dynamic_range()
        ax[1].scatter(1E-6*stat_KA['mode'].values[0],1E-6*stat_KI['mode'].values[0],
                    color = colors[i[0]] , edgecolor='white')

divider = make_axes_locatable(ax[1])
cax = divider.append_axes("right", size="5%", pad=0.2)
cax.set_aspect(12)

v = np.arange(0, 1.2, 0.2)
clb = plt.colorbar(sc1, cax=cax, ticks=v )
clb.ax.set_title('dynamic\nrange',loc='center',position=(1.6, 1.05), fontsize=8)

ax[1].set_yscale('log')
ax[1].set_xscale('log')
ax[1].set_xlabel('$K_A$ (M)', fontsize=8)
ax[1].set_ylabel('$K_I$ (M)', fontsize=8)
ax[1].tick_params(labelsize=8)

# for plotting the legend (removing replicate labels)
handles, labels = ax[0].get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
ax[0].legend(by_label.values(), by_label.keys(),
             loc='upper left', fontsize=8)

plt.tight_layout()

fig.savefig('fig4.pdf', bbox_inches='tight')
