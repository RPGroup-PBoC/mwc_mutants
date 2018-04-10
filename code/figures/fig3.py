#%%
# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import pandas as pd
sys.path.insert(0, '../../')
import mut.thermo
import mut.bayes
import mut.stats
import mut.viz
colors = mut.viz.pub_style()
FIG_NO = 3

# Load the necessary data sets.
# data_dir = '../data'
data_dir = 'data/'
titration_data = pd.read_csv("{}merged_dna_O2_data.csv".format(data_dir))
epR_fit_params = pd.read_csv("{}epR_fit_DNA_O2_stats.csv".format(data_dir))
global_fit_stats = pd.read_csv('{}global_fit_DNA_O2_stats.csv'.format(data_dir))
ka_lit = 139E-6
ki_lit = 0.53E-6
wt_epR = -13.9
wt_R = 260
c_range = np.logspace(-8, -2, 500)
rep_range = np.logspace(0, 4, 500)
wt_titration = mut.thermo.SimpleRepression(wt_R, wt_epR, effector_conc=c_range, ep_ai=4.5, ka=ka_lit, ki=ki_lit,
n_sites = 2)
wt_leakiness = mut.thermo.SimpleRepression(rep_range, wt_epR, effector_conc=0, ep_ai=4.5, ka=ka_lit, ki=ki_lit,
n_sites=2)
wt_fc = wt_titration.fold_change()
wt_leak = wt_leakiness.fold_change()
color_dict = {1220:colors['blue'], 260:colors['red'], 124:colors['yellow'], 60:colors['green']}


#  Compute the mean and sem for the titration data
def compute_mean_sem(df):
    mean_fc = df['fold_change'].mean()
    sem_fc = df['fold_change'].std() / np.sqrt(len(df))
    return mean_fc, sem_fc 

grouped = titration_data.groupby(['mutant', 'repressors', 'IPTGuM'])

# %%
import imp
imp.reload(mut.thermo)
# # Set up the figure axis.
fig = plt.figure(figsize=(6, 6))
gs = gridspec.GridSpec(10, 12)

# Leakiness vals
ax0 = fig.add_subplot(gs[0:3, 0:4])
ax1 = fig.add_subplot(gs[0:3, 4:8])
ax2 = fig.add_subplot(gs[0:3, 8:])

# Titrations
ax3 = fig.add_subplot(gs[3:6, 0:4])
ax4 = fig.add_subplot(gs[3:6, 4:8])
ax5 = fig.add_subplot(gs[3:6, 8:])

# Ka/Ki Collapse
ax6 = fig.add_subplot(gs[6:, 0:4])
ax8 = fig.add_subplot(gs[6:, 4:])
ax = [ax0, ax1, ax2, ax3, ax4, ax5, ax6]

ax0.axis('off')
ax1.set_xlabel('repressors per cell')
ax2.set_ylabel('binding energy\n [ $k_BT$ ]')

ax8.set_ylabel('fold-change')
ax8.set_xlabel('Bohr parameter [ $k_BT$ ]')

# Format the axes as necessary
ax1.set_xscale('log')
ax1.set_yscale('log')
for a in ax[3:6]:
    a.set_ylim([-0.1, 1.2])
    a.set_xlim([1E-8, 1E-2])
    a.set_xscale('log')
    a.set_xlabel('IPTG [M]')
    a.set_ylabel('fold-change')
mut_axes = {'Q21M':ax[3], 'Q21A':ax[4], 'Y20I':ax[5]}

# Plot the WT controls
_ = ax1.plot(rep_range, wt_leak, 'k', lw=0.5, label='WT, 260')
for a in ax[3:6]:
    _ = a.plot(c_range, wt_fc, 'k', lw=0.5, label='WT, 260')
_ = ax3.hlines(wt_epR, 0, 4, 'k', lw=0.5)

# Get the mode and hpd of the fit. 
vals = {}
for i, m in enumerate(titration_data['mutant'].unique()): 
    if m != 'wt':
        vals[m] = epR_fit_params[epR_fit_params['parameter'] == 'epR_{}'.format(m)].values[0][1:]

# Plot the leakiness for each mutant
for i, m in enumerate(titration_data['mutant'].unique()):
    if m != 'wt':
        R_mesh,  epR_mesh = np.meshgrid(rep_range, list(vals[m]))
        arch = mut.thermo.SimpleRepression(R_mesh, epR_mesh, ka=ka_lit, ki=ki_lit, ep_ai=4.5, 
                                       effector_conc=0, n_sites=2)
        leak = arch.leakiness()

        R_range = titration_data['repressors'].unique()
        R_mesh,  c_mesh, epR_mesh = np.meshgrid(R_range, c_range, list(vals[m]))
        arch = mut.thermo.SimpleRepression(R_mesh, epR_mesh, ka=ka_lit, ki=ki_lit, ep_ai=4.5, 
                                       effector_conc=c_mesh, n_sites=2)
        fc = arch.fold_change()
        for k, r in enumerate(R_range): 
            if m=='Q21M':
                label = int(r)
            else:
                label = '__nolegend__'
            _ = mut_axes[m].plot(c_range, fc[:, k, 0], color=color_dict[int(r)], lw=1, label=label)
            _ = mut_axes[m].fill_between(c_range, fc[:, k, 1], fc[:, k, 2], color=color_dict[int(r)], alpha=0.5)
            _ = mut_axes[m].set_title(m, backgroundcolor=colors['pale_yellow'], y=1.01, fontsize=8)

        _ = ax1.plot(rep_range, leak[0, :], color=leak_dict[m], lw=1)
        _ = ax1.fill_between(rep_range, leak[1, :], leak[2, :], color=leak_dict[m], alpha=0.5)

# Plot the theoretical bohr. 
f_range = np.linspace(-10, 10)
bohr_theo = (1 + np.exp(-f_range))**-1
_ = ax8.plot(f_range, bohr_theo, color=colors['blue'], linewidth=2)

#Plot the data
viridis = sns.color_palette('viridis', n_colors=3)
leak_dict = {'Q21M': viridis[0], 'Q21A': viridis[1], 
            'Y20I': viridis[2], 'wt':'slategrey'} 
# Leakiness
for g, d in grouped:
    mean_fc, sem_fc = compute_mean_sem(d)
    if g[-1] == 0:
        if g[0] == 'wt':
            pass            
        else:
            _ = ax[1].errorbar(g[1], mean_fc, sem_fc, color=leak_dict[g[0]], fmt='o',
                             lw=1, markersize=2)


# Induction curves
for g, d in grouped:
    mean_fc, sem_fc = compute_mean_sem(d)
    if g[0] == 'wt':    
        pass
    else: 
        _ = mut_axes[g[0]].errorbar(g[-1]/1e6, mean_fc, sem_fc,
                color= color_dict[g[1]], fmt='o', lw=1, markersize=2)


# Values for epR fit.
arrangement = {'Q21M': 1, 'Q21A':2, 'Y20I': 3}
modes = {}
for i, m in enumerate(titration_data['mutant'].unique()): 
    param = epR_fit_params[epR_fit_params['parameter']=='epR_{}'.format(m)]
    modes[m] = param['mode'].values[0]
    if m == 'wt':
        pass
    else:
        _ = ax2.plot(arrangement[m], param['mode'], 'o', color=leak_dict[m], ms=2)
        _ = ax2.vlines(arrangement[m], param['hpd_min'], param['hpd_max'], lw=1, color=leak_dict[m])

# Plot the Hernan determined value for O2
_ = ax2.hlines(-13.9, 0, 4, color='slategray', lw=1)

# Plot Ka/Ki
for i, m in enumerate(titration_data['mutant'].unique()):
    if m != 'wt':
        ka_param = global_fit_stats[global_fit_stats['parameter']=='ka_{}'.format(m)]
        ki_param = global_fit_stats[global_fit_stats['parameter']=='ki_{}'.format(m)]
        _ = ax6.plot(arrangement[m], ka_param['mode'].values[0] / ki_param['mode'].values[0], 'o', color=leak_dict[m])
        _ = ax6.vlines(arrangement[m], ka_param['hpd_max'].values[0] / ki_param['hpd_min'].values[0],
                        ka_param['hpd_min'].values[0] / ki_param['hpd_max'].values[0], lw=1, color=leak_dict[m])

ax6.hlines(ka_lit / ki_lit, 0, 5, 'k', lw=0.5)
ax6.set_xlim([0, 4])
ax6.set_xticks([0, 1, 2, 3, 4])
ax6.set_xticklabels(['', 'Q21M', 'Q21A', 'Y20I', ''], rotation=45)
ax6.xaxis.grid(False)
ax6.set_yscale('log')
# Data collapse
# Set up  the architecture.
for g, d in grouped:
    if g[0] == 'wt':
        pass
    else:
        mean_fc, sem_fc = compute_mean_sem(d)
        arch = mut.thermo.SimpleRepression(R=g[1], ep_r=modes[g[0]], effector_conc=g[-1]/1E6,
        ka=ka_lit, ki=ki_lit, ep_ai=4.5, n_sites=2)
        bohr_param = arch.bohr_parameter()

        # Plot
        _ = ax8.errorbar(bohr_param, mean_fc, sem_fc, fmt='o', color=leak_dict[g[0]], markersize=2,
        lw=1)


ax2.xaxis.grid(False)
ax2.set_xticks([1, 2, 3])
ax2.set_yticks([-16, -14, -12, -10, -8])
ax2.set_yticklabels([-16, -14, -12, -10, -8])

ax2.set_xlim([0, 4])
ax2.set_xticklabels(arrangement.keys(), rotation=45)
for a in ax[3:6]:
     a.set_xticks([1E-8, 1E-5, 1E-2])

leg =  ax3.legend(loc='upper left', fontsize=6, title='rep. / cell', handlelength=1)
plt.setp(leg.get_title(), fontsize=6)
plt.tight_layout()
# plt.subplots_adjust(hspace=2, wspace=80)
plt.savefig('fig.pdf')

