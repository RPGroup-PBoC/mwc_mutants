# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.bayes
import mut.thermo
import mut.viz
import seaborn as sns
pboc = mut.viz.color_selector('pboc')
color_palette = sns.color_palette('deep', n_colors=18)
mut.viz.plotting_style()

constants = mut.thermo.load_constants()
data = pd.read_csv('../../data/csv/RazoMejia2018_data.csv')
data = data[data['repressors'] > 0]
data['repressors'] *= 2
data['mutant'] = 'wt'
data['class'] = 'wt'
data.rename(columns={'fold_change_A':'fold_change', 
                    'IPTG_uM':'IPTGuM'}, inplace=True)

rep_colors = {22:pboc['red'], 60:pboc['blue'], 124:pboc['green'],
             260:pboc['purple'], 1220:pboc['dark_green'], 1740:pboc['dark_brown']}

operator_glyphs = {'O1': 's', 'O2': 'o', 'O3': '^'}

# Define the reference states.
c0 = 50; # in µM
R0 = constants['RBS1027']
epRA0 = constants['O2']
data
R_ref = mut.thermo.SimpleRepression(R=R0, ep_r=data['binding_energy'],
                                   ka=constants['Ka'], ki=constants['Ki'],
                                   n_sites=constants['n_sites'],
                                   effector_conc=data['IPTGuM'],
                                   ep_ai=constants['ep_AI']).bohr_parameter()
ep_ref = mut.thermo.SimpleRepression(R=data['repressors'], ep_r=epRA0,
                                   ka=constants['Ka'], ki=constants['Ki'],
                                   n_sites=constants['n_sites'],
                                   effector_conc=data['IPTGuM'],
                                   ep_ai=constants['ep_AI']).bohr_parameter()
c_ref = mut.thermo.SimpleRepression(R=data['repressors'], ep_r=data['binding_energy'],
                                   ka=constants['Ka'], ki=constants['Ki'],
                                   n_sites=constants['n_sites'],
                                   effector_conc=c0,
                                   ep_ai=constants['ep_AI']).bohr_parameter()
# Compute the statistics. 
data['ref_bohr'] = R_ref
R_stats = mut.bayes.infer_empirical_bohr(data, '../stan/empirical_F.stan')
data['ref_bohr'] = ep_ref
ep_stats = mut.bayes.infer_empirical_bohr(data, '../stan/empirical_F.stan')
data['ref_bohr'] = c_ref
c_stats = mut.bayes.infer_empirical_bohr(data, '../stan/empirical_F.stan')


# Define the parameter ranges for the collapse curves 
r_range = np.logspace(0, 4, 200)
ep_range = np.linspace(-17, -2, 200)
c_range = np.logspace(-4, 4, 200)
pact_ref = mut.thermo.MWC(ka=constants['Ka'], ki=constants['Ki'],
                         n_sites=constants['n_sites'], ep_ai=constants['ep_AI'],
                         effector_conc=c0).pact()
pact_range = mut.thermo.MWC(ka=constants['Ka'], ki=constants['Ki'],
                         n_sites=constants['n_sites'], ep_ai=constants['ep_AI'],
                         effector_conc=c_range).pact()

# Compute the analytical collapse curves
ep_collapse = epRA0 - ep_range
r_collapse = -np.log(R0/r_range)
pact_collapse = -np.log(pact_ref / pact_range)


# Instantiate the figure
fig, ax = plt.subplots(1, 3, figsize=(6, 1.5))

for g, d in c_stats.groupby(['repressors', 'operator']):
    _d = d[d['parameter']=='delta_bohr_corrected']
    pact = mut.thermo.MWC(ka=constants['Ka'], ki=constants['Ki'],
                                      ep_ai=constants['ep_AI'], n_sites=constants['n_sites'],
                                      effector_conc=_d['IPTGuM']).pact()
    _ = ax[0].plot(pact_ref/pact, _d['median'], marker='o',
                  color='slategray', ms=3, linestyle='none', alpha=0.5,
                  )
    
for g, d in R_stats.groupby(['repressors', 'operator']):
    _d = d[d['parameter']=='delta_bohr_corrected']
    _ = ax[1].plot(R0/_d['repressors'], _d['median'], marker='o',
                  color='slategray', ms=3, linestyle='none', 
                  alpha=0.5)

for g, d in ep_stats.groupby(['repressors', 'operator']):
    _d = d[d['parameter']=='delta_bohr_corrected']
    _ = ax[2].plot(epRA0 - constants[g[1]] * np.ones(len(_d)), _d['median'], marker='o',
                  color='slategray', ms=3, linestyle='none', alpha=0.5,
                  markerfacecolor=None)
    
c_grouped = c_stats[c_stats['parameter']==\
                    'delta_bohr_corrected'].groupby(['IPTGuM']).agg(
                    ('mean', 'std')).reset_index() 
R_grouped = R_stats[R_stats['parameter']==\
                    'delta_bohr_corrected'].groupby(['repressors']).agg(
                    ('mean', 'std')).reset_index()
ep_grouped = ep_stats[ep_stats['parameter']==\
                    'delta_bohr_corrected'].groupby(['operator']).agg(
                    ('mean', 'std')).reset_index()
ep_binding_energy = np.array([constants[op] for op in ep_grouped['operator'].values])
pact = mut.thermo.MWC(ka=constants['Ka'], ki=constants['Ki'], n_sites=constants['n_sites'],
                     ep_ai=constants['ep_AI'],  effector_conc=c_grouped['IPTGuM']).pact()
_ = ax[0].errorbar(pact_ref/pact, c_grouped['median']['mean'], c_grouped['median']['std'],
                   linestyle='none', color='black', fmt='.', lw=1, zorder=1000,
                  markerfacecolor='w', capsize=2)
_ = ax[1].errorbar(R0/R_grouped['repressors'], R_grouped['median']['mean'], R_grouped['median']['std'],
                   linestyle='none', color='black', fmt='.', lw=1, zorder=1000,
                  markerfacecolor='w', capsize=2)
_ = ax[2].errorbar(epRA0 - ep_binding_energy, ep_grouped['median']['mean'], ep_grouped['median']['std'],
                   linestyle='none', color='black', fmt='.', lw=1, zorder=1000,
                  markerfacecolor='w', capsize=2)  
ax[0].plot(pact_ref/pact_range, pact_collapse, color=pboc['blue'], lw=2)
ax[1].plot(R0/r_range, r_collapse, pboc['red'], lw=2)
ax[2].plot(epRA0 - ep_range, ep_collapse, pboc['purple'], lw=2)

# Add appropriate labels. 
ax[0].set_xlabel('$p_\mathrm{act}^{(\mathrm{ref})}(c)) / p_\mathrm{act}^*(c)$', fontsize=8)
ax[1].set_xlabel('$R^{(\mathrm{ref})} / R^*$', fontsize=8)
ax[2].set_xlabel(r'$\Delta\varepsilon_{RA}^{(\mathrm{ref})} - \Delta\varepsilon_{RA}^*$ [$k_BT$]', fontsize=8)

for a in ax:
    a.set_ylabel('$F^{(\mathrm{ref})} - F^*$ [$k_BT$]', fontsize=8)
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.grid(False)
    
# Set the various limits    
ax[0].set_xlim([1E-2, 1E1])
ax[1].set_xlim([R0/1, R0/1E4])
ax[2].set_xlim((-7,3))
ax[2].set_ylim([-10, 5])
ax[0].set_xscale('log')
ax[1].set_xscale('log')
# plt.tight_layout()
plt.savefig('../../figures/fig1_RazoMejia2018_deltaF.svg', bbox_inches='tight')


