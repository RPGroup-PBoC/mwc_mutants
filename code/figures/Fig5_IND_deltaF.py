# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.thermo
import mut.stats
import mut.viz
import seaborn as sns
constants = mut.thermo.load_constants()
_colors = sns.color_palette('deep')
op_colors = {'O1': _colors[0], 'O2': _colors[1], 'O3':_colors[2]}
pboc = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load the data, 
deltaF_data = pd.read_csv('../../data/csv/empirical_F_statistics.csv')
data = deltaF_data[deltaF_data['class']=='IND']

# Load the statistics
kaki_only_stats = pd.read_csv('../../data/csv/KaKi_only_summary.csv')
kaki_only_stats = kaki_only_stats[kaki_only_stats['operator']=='O2']
kaki_epai_stats = pd.read_csv('../../data/csv/KaKi_epAI_summary.csv')
kaki_epai_stats = kaki_epai_stats[kaki_epai_stats['operator']=='O2']
kaki_epai_samps = pd.read_csv('../../data/csv/KaKi_epAI_samples.csv')
kaki_epai_samps = kaki_epai_samps[kaki_epai_samps['operator']=='O2']

# Experimental constants
c_range = np.logspace(-8, 4, 100)

# ##############################################################################
# FIGURE INSTANTIATION AND FORMATTING
# ##############################################################################
fig, ax = plt.subplots(4, 3, figsize=(3.42, 4), sharey=True)

op_ind = {'O1':0, 'O2':1, 'O3':2}
mut_ind = {'F164T': 3, 'Q294V': 2, 'Q294K':1, 'Q294R': 0}
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    # a.set_xlim([-8, 8])

for o, ind in op_ind.items():
    ax[0, ind].set_title(o, fontsize=7, backgroundcolor=pboc['pale_yellow'], y=1.08)


# ##############################################################################
# DELTA F DATA
# ##############################################################################
for g, d in data.groupby(['mutant', 'operator']):
    param = d[d['parameter']=='delta_bohr_corrected'] 
    ref = mut.thermo.SimpleRepression(R=260, ep_r=constants[g[1]],
                                      ka=constants['Ka'], ki=constants['Ki'],
                                      ep_ai=constants['ep_AI'],
                                      effector_conc=param['IPTGuM']).bohr_parameter()
    _ax = ax[mut_ind[g[0]], op_ind[g[1]]]

    # Determine conditional coloring
    if g[1] == 'O2':
        face='w'
    else:
        face = op_colors[g[1]]

    _ax.plot(ref, param['median'], 'o', markerfacecolor=face, color=op_colors[g[1]],
            ms=2)
    _ax.vlines(ref, param['hpd_min'], param['hpd_max'], color=op_colors[g[1]], lw=1)
    
# ##############################################################################
# DELTA F CURVES
# ##############################################################################
for i, m in enumerate(mut_ind.keys()):
    for j, o in enumerate(op_ind.keys()):
        # Compute the reference free energy
        ref = mut.thermo.SimpleRepression(R=260, ep_r=constants[o],
                                          ka=constants['Ka'], ki=constants['Ki'], 
                                         ep_ai=constants['ep_AI'],
                                          effector_conc=c_range).bohr_parameter()

        # Get the statistics for the  Ka Ki only fits.
        _stats = kaki_only_stats[kaki_only_stats['mutant']==m]
        ka_median = _stats[_stats['parameter']=='Ka']['median'].values[0]
        ki_median = _stats[_stats['parameter']=='Ki']['median'].values[0]
        fit1 = mut.thermo.SimpleRepression(R=260, ep_r=constants[o],
                                            ka=ka_median, ki=ki_median, 
                                              ep_ai=constants['ep_AI'],
                                              effector_conc=c_range).bohr_parameter()

        _stats = kaki_epai_stats[kaki_epai_stats['mutant']==m]
        _samps = kaki_epai_samps[kaki_epai_samps['mutant']==m]
        ka_median = _stats[_stats['parameter']=='Ka']['median'].values[0]
        ki_median = _stats[_stats['parameter']=='Ki']['median'].values[0]
        epai_median = _stats[_stats['parameter']=='ep_AI']['median'].values[0]
        fit = mut.thermo.SimpleRepression(R=260, ep_r=constants[o],
                             ka=ka_median, ki=ki_median, 
                             ep_ai=epai_median,
                             effector_conc=c_range).bohr_parameter()
        cred_region = np.zeros((2, len(c_range)))
        for k, c in enumerate(c_range):
            arch = mut.thermo.SimpleRepression(R=260, ep_r=constants[o],
                                                ka=_samps['Ka'], ki=_samps['Ki'],
                                                ep_ai=_samps['ep_AI'],
                                                effector_conc=c).bohr_parameter()
            delF = ref[k] - arch
            cred_region[:, k] = mut.stats.compute_hpd(delF, 0.95)

        # ax[mut_ind[m], op_ind[o]].plot(ref, ref - fit, '--', color=op_colors[o])
        ax[mut_ind[m], op_ind[o]].fill_between(ref, cred_region[0, :],
                                               cred_region[1, :], color=op_colors[o],
                                               alpha=0.4)

    

plt.subplots_adjust(hspace=0.1, wspace=0.1)
plt.savefig('/Users/gchure/Desktop/ind_muts.pdf')



