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


# Define the parameter ranges for the collapse curves 
r_range = np.logspace(0, 4, 200)
ep_range = np.linspace(-17, -2, 200)
c_range = np.logspace(-2, 4, 200)
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

# Compute the statistics. 
data['ref_bohr'] = R_ref
R_stats = mut.bayes.infer_empirical_bohr(data, '../stan/empirical_F.stan')
data['ref_bohr'] = ep_ref
ep_stats = mut.bayes.infer_empirical_bohr(data, '../stan/empirical_F.stan')
data['ref_bohr'] = c_ref
c_stats = mut.bayes.infer_empirical_bohr(data, '../stan/empirical_F.stan')


# Instantiate the figure
fig, ax = plt.subplots(1, 3, figsize=(7, 2))

iter = 0
for g, d in c_stats.groupby(['repressors', 'operator']):
    _d = d[d['parameter']=='delta_bohr_corrected']
    _ = ax[0].plot(_d['IPTGuM']/c0, _d['median'], marker='o',
                  color='slategray', ms=2, linestyle='none', alpha=0.5,
                  markerfacecolor='none')
    iter += 1


iter = 0
for g, d in R_stats.groupby(['repressors', 'operator']):
    _d = d[d['parameter']=='delta_bohr_corrected']
    _ = ax[1].plot(_d['repressors']/R0, _d['median'], marker='o',
                  color='slategray', ms=2, linestyle='none', 
                  alpha=0.5)

iter = 0
for g, d in ep_stats.groupby(['repressors', 'operator']):
    _d = d[d['parameter']=='delta_bohr_corrected']
    _ = ax[2].plot(epRA0 - constants[g[1]] * np.ones(len(_d)), _d['median'], marker='o',
                  color='slategray', ms=2, linestyle='none', alpha=0.5,
                  markerfacecolor=None)
    iter += 1
    
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
_ = ax[0].errorbar(c_grouped['IPTGuM']/c0, c_grouped['median']['mean'], c_grouped['median']['std'],
                   linestyle='none', color='black', fmt='.', lw=1, zorder=1000,
                  markerfacecolor='w', capsize=2)
_ = ax[1].errorbar(R_grouped['repressors']/R0, R_grouped['median']['mean'], R_grouped['median']['std'],
                   linestyle='none', color='black', fmt='.', lw=1, zorder=1000,
                  markerfacecolor='w', capsize=2)
_ = ax[2].errorbar(epRA0 - ep_binding_energy, ep_grouped['median']['mean'], ep_grouped['median']['std'],
                   linestyle='none', color='black', fmt='.', lw=1, zorder=1000,
                  markerfacecolor='w', capsize=2)  
ax[0].plot(c_range/c0, pact_collapse, color=pboc['blue'], lw=2)
ax[1].plot(r_range/R0, r_collapse, pboc['red'], lw=2)
ax[2].plot(epRA0 - ep_range, ep_collapse, pboc['purple'], lw=2)


# Set the various limits    
ax[0].set_xscale('log')
ax[1].set_xscale('log')
# ax[2].set_xlim([-1, 1E-2])
plt.savefig('/Users/gchure/Desktop/test.pdf')


