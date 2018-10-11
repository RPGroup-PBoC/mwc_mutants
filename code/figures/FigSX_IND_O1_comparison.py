# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
sys.path.insert(0, '../../')
import mut.thermo
import mut.viz
import mut.stats
pboc = mut.viz.color_selector('pboc')
colors = mut.viz.color_selector('mut')
constants = mut.thermo.load_constants()
mut.viz.plotting_style()

# Load and restrict the summarized data
data = pd.read_csv('../../data/csv/summarized_data.csv')
data = data[((data['mutant']=='Q294K') | (data['mutant']=='Q294R')) & (data['operator'] != 'O3')]

# Load the mcm samples.
O2_samples = pd.read_csv('../../data/csv/Fig3_O2_KaKi_epAI_samples.csv')
O2_modes = O2_samples.iloc[np.argmax(O2_samples['logp'].values)]
O1_samples = pd.read_csv('../../data/csv/Fig3_O1_KaKi_epAI_samples.csv')
O1_modes = O1_samples.iloc[np.argmax(O1_samples['logp'].values)]
modes = {'O1':O1_modes, 'O2':O2_modes}
samples = {'O1':O1_samples, 'O2':O2_samples}

# Define parameter ranges. 
c_range = np.logspace(-2, 4, 200)

# Set up the figure canvas. 
fig, axes = plt.subplots(2, 3, figsize=(6, 4.7))
ax = [axes[0,0], axes[0, 1], axes[1, 0], axes[1, 1]]

# Properly format each axis. 
for a in axes.ravel():
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
for i, a in enumerate(ax): 
    a.set_xscale('log')
    a.set_xlim([1E-8, 1E-2])
    a.set_xticks([1E-7, 1E-5, 1E-3])
    a.set_ylim([-0.05, 1.2]) 
    if (i == 0) | (i == 3):
        a.set_facecolor('#e4e7ec')

# Add conditional formatting                
ax[0].set_xticklabels([])
ax[1].set_xticklabels([])
ax[1].set_yticklabels([])
ax[3].set_yticklabels([])
axes[0, 2].set_xlim([-3.75, -1.5])
axes[1, 2].set_xlim([-3.75, -1.5])
axes[0, 2].set_xticklabels([])
axes[0, 2].set_yticks([0])
axes[1, 2].set_yticks([0])
axes[0, 2].set_yticklabels([])
axes[1, 2].set_yticklabels([])


# Add axis labels
ax[0].set_ylabel('fold-change', fontsize=8)
ax[2].set_ylabel('fold-change', fontsize=8)
ax[2].set_xlabel('IPTG [M]', fontsize=8)
ax[3].set_xlabel('IPTG [M]', fontsize=8)
axes[0, 2].set_title('posterior distributions\n for $\Delta\varepsilon_{AI}$', fontsize=8,
                    backgroundcolor=pboc['pale_yellow'], y=1.04)

# Add titles and labels
ax[0].set_title(r'O1 | $\Delta\varepsilon_{RA}=%s\, k_BT$' %constants['O1'], fontsize=8,
               backgroundcolor=pboc['pale_yellow'], y=1.04)
ax[1].set_title(r'O2 | $\Delta\varepsilon_{RA}=%s\, k_BT$' %constants['O2'], fontsize=8,
               backgroundcolor=pboc['pale_yellow'], y=1.04)
ax[0].text(-0.49, 0.75, r'O1 | $\Delta\varepsilon_{RA}=%s\, k_BT$' %constants['O1'], fontsize=8,
          rotation='vertical', backgroundcolor=pboc['pale_yellow'], transform=ax[0].transAxes)
ax[2].text(-0.49, 0.75, r'O2 | $\Delta\varepsilon_{RA}=%s\, k_BT$' %constants['O2'], fontsize=8,
          rotation='vertical', backgroundcolor=pboc['pale_yellow'], transform=ax[2].transAxes)
axes[1, 2].set_xlabel(r'$\Delta\varepsilon_{AI}$ [$k_BT$]', fontsize=8)

# Plot the summarized data. 
op_ax = {'O1':[ax[0], ax[2]], 'O2': [ax[3], ax[1]]}
for g, d in data.groupby(['operator', 'mutant']): 
    # Plot the self data. 
    _ = op_ax[g[0]][0].errorbar(d['IPTGuM'] / 1E6, d['mean'], d['sem'],
                         markerfacecolor='w', markersize=4, markeredgecolor=colors[g[1]],
                        lw=1, linestyle='none', fmt='o', color=colors[g[1]], label=g[1])

    # Plot the comparison data
    _ = op_ax[g[0]][1].errorbar(d['IPTGuM'] / 1E6, d['mean'], d['sem'],
                               color=colors[g[1]], markersize=4, label=g[1], linestyle='none',
                               fmt='o', lw=1)
    

pre_op = {'O1':[constants['O1'], constants['O2']], 'O2':[constants['O1'], constants['O2']]}   
fit_ax = {'O1':[ax[0], ax[1]], 'O2':[ax[2], ax[3]]}
op_colors = {'O1':pboc['red'], 'O2':pboc['blue']}
for m in data['mutant'].unique(): 
    for o, e in pre_op.items():
        _modes = modes[o]
        _samples = samples[o]
        
        # Plot the posterior distribution for ep_AI for the two mutants.
        if o == 'O1':
            _i = 0
        else:
            _i = 1 
        _ = axes[_i, 2].hist(_samples[f'ep_AI.{m}'], histtype='stepfilled', color=colors[m],
                           alpha=0.9, edgecolor=colors[m], lw=1, bins=100, label=o)
    
        for i in range(2):
            # Find the parameters and plot
            arch = mut.thermo.SimpleRepression(R=d['repressors'].unique()[0], ep_r=e[i],
                                          effector_conc=c_range, ep_ai=_modes[f'ep_AI.{m}'],
                                          ka=_modes[f'Ka.{m}'], ki=_modes[f'Ki.{m}'],
                                          n_ns=constants['Nns'], n_sites=constants['n_sites']).fold_change()
        
            # Compute the credible regions for the self fit. 
            cred_regions = np.zeros([2, len(c_range)])
            for j, c in enumerate(c_range):
                _arch = mut.thermo.SimpleRepression(R=d['repressors'].unique()[0], ep_r=e[i],
                                              effector_conc=c, ep_ai=_samples[f'ep_AI.{m}'],
                                              ka=_samples[f'Ka.{m}'], ki=_samples[f'Ki.{m}'],
                                              n_ns=constants['Nns'], n_sites=constants['n_sites']).fold_change()
                cred_regions[:, j] = mut.stats.compute_hpd(_arch, mass_frac=0.95)
          
            _ = fit_ax[o][i].plot(c_range/1E6, arch, lw=1, color=colors[m], label='__nolegend__')
            _ = fit_ax[o][i].fill_between(c_range/1E6, cred_regions[0, :], cred_regions[1, :], color=colors[m],
                                            alpha=0.5, label='__nolegend__')




fig.text(0.3, 1.02, 'comparison strain', backgroundcolor='#E3DCD0', fontsize=10)   
fig.text(-0.05, 0.54, 'fit strain', backgroundcolor='#E3DCD0', fontsize=10, rotation='vertical') 
# plt.subplots_adjust(wspace=0.04, hspace=0.04)
plt.tight_layout()
ax[1].legend(fontsize=8)
plt.savefig('FigSX_IND_pairwise_prediction.pdf', bbox_inches='tight')