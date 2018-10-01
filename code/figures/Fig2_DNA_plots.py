# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.thermo
import mut.viz
pboc = mut.viz.color_selector('pboc')
colors = mut.viz.color_selector('mut')
constants = mut.thermo.load_constants()
mut.viz.plotting_style()

# Load the DNA binding domain mutants. 
data = pd.read_csv('../../data/csv/summarized_data.csv')
data = data[((data['class']=='DNA') | (data['class']=='WT')) & (data['operator']=='O2')].copy()

# Load the fit parameters.
epRA_stats = pd.read_csv('../../data/csv/Fig2_O2_DNA_binding_energy_stats.csv')

# Define the plot constants. 
c_range = np.logspace(-8, -2, 200)
rep_range = np.logspace(0, 4, 200)

# Define the glyphs for plotting based on repressor copy number. 
_glyphs = ['o', 's', 'D', '^']
glyphs = {i:j for i, j in zip(data['repressors'].unique(), _glyphs)}

# Define functions for computing the dynamic range as a function of R / kdna
def dynamic_range(r_kdna, constants):
    allo = mut.thermo.MWC(ep_ai=constants['ep_AI'], ka=ka, ki=ki, n_sites=constants['n_sites'], effector_conc=0)
    sat = (1 + allo.saturation() * r_kdna)**-1
    leak = (1 + allo.leakiness() * r_kdna)**-1
    return sat - leak

# Instantiate the figure
fig, ax = plt.subplots(2, 2, figsize=(5, 4.5))

# Format the axes and add appropriate labels.
ax[0, 0].axis('off')
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_ylim([-0.05, 1.2])
ax[0, 1].set_xscale('log')
ax[1, 0].set_xscale('log')
ax[0, 1].set_xlabel('IPTG [M]', fontsize=8)
ax[0, 1].set_ylabel('fold-change', fontsize=8)
ax[1, 0].set_ylabel('dynamic range', fontsize=8)
ax[1, 0].set_xlabel(r'$\frac{R}{N_{NS}}e^{-\beta\Delta\varepsilon_{RA}}$', fontsize=8)
ax[1, 1].set_ylabel('fold-change', fontsize=8)
ax[1, 1].set_xlabel('Bohr parameter [$k_BT$]', fontsize=8)

# Restrict the axes. 
ax[0, 1].set_xlim([1E-8, 1E-2])
ax[1, 1].set_xlim([-9, 9])
ax[1, 0].set_xlim([1E-2, 1E5])
#Add panel labels.
labels = ['A', 'B', 'C', 'D']
for i, a in enumerate(ax.ravel()):
    a.text(-0.25, 1.05, '({})'.format(labels[i]), fontsize=8,
          transform=a.transAxes)

# Compute the predicted dynamic range
r_kdna_range = np.logspace(-2, 5, 200)
dyn_rng_pred = dynamic_range(r_kdna_range, constants)
ax[1, 0].plot(r_kdna_range, dyn_rng_pred, color='#3c3c3c', lw=1)

# Compute the data collapse curve. 
bohr_range = np.linspace(-9, 9, 200)
collapse = (1 + np.exp(-bohr_range))**-1
ax[1, 1].plot(bohr_range, collapse, lw=1, color='#3c3c3c', zorder=1000)
    
# Add the glyphs
for g, d in data.groupby(['mutant', 'repressors']): 
    # Isolate the statistics
    if g[0].lower() != 'wt':
        epRA = epRA_stats[epRA_stats['parameter']=='ep_RA.{}.260'.format(g[0])].values[0][1:].astype(float)
    else:
        epRA = np.array([constants['O2']] * 3)
    ka = constants['Ka'] / 1E6
    ki = constants['Ki'] / 1E6
        
    # Compute the observed dynamic range and plot. 
    dyn_rng = d[d['IPTGuM']==5000]['mean'].values - d[d['IPTGuM']==0]['mean'].values  
    r_kdna = (g[1] / constants['Nns']) * np.exp(-epRA)
    ax[1, 0].plot(r_kdna[0], dyn_rng, marker=glyphs[g[1]], ms=3, color=colors[g[0].upper()])
    ax[1, 0].hlines(dyn_rng, r_kdna[1], r_kdna[2], lw=1, color=colors[g[0].upper()])
    
    # Compute the bohr parameter.
    c_mesh, epRA_mesh = np.meshgrid(d['IPTGuM'], epRA)
    bohr = mut.thermo.SimpleRepression(R=g[1], ep_r=epRA_mesh, effector_conc=c_mesh,
                                      ka=constants['Ka'], ki=constants['Ki'], 
                                      ep_ai=constants['ep_AI'], n_sites=constants['n_sites'],
                                      n_ns=constants['Nns']).bohr_parameter()
    ax[1, 1].errorbar(bohr[0, :], d['mean'], d['sem'], color=colors[g[0].upper()], fmt='o', 
                          ms=2, linestyle='none', lw=1, label='__nolegend__')
    ax[1, 1].hlines(d['mean'], bohr[1,:], bohr[2, :], color=colors[g[0].upper()], lw=1)
    
    # Plot the induction profiles if repressors = 260. 
    if g[1] == 260: 
        ax[0, 1].errorbar(d['IPTGuM'] / 1E6, d['mean'], d['sem'], 
                         fmt='o', linestyle='none', lw=1, ms=3,
                         color=colors[g[0].upper()], label=g[0].upper())
        
        # Compute the fits and plot. 
        c_mesh, epRA_mesh = np.meshgrid(c_range, epRA)
        fits = mut.thermo.SimpleRepression(R=g[1], ep_r=epRA_mesh, ka=ka, ki=ki,
                                          ep_ai=constants['ep_AI'], effector_conc=c_mesh, 
                                          n_sites=constants['n_sites'], 
                                          n_ns=constants['Nns']).fold_change()
        ax[0, 1].plot(c_range, fits[0, :], color=colors[g[0].upper()], lw=1, 
                     label='__nolegend__')
        ax[0, 1].fill_between(c_range, fits[1, :], fits[2, :], alpha=0.3, color=colors[g[0].upper()])

# Add legends
for r, g in glyphs.items():
    ax[1, 0].plot([], [], linestyle='none', color='k', alpha=0.3, marker=g, label=int(r),
                 ms=3)
_leg = ax[1, 0].legend(bbox_to_anchor=(0.7, 0.98), fontsize=7.5, handletextpad=0.07,
                      title='rep. / cell')
_leg.get_title().set_fontsize(7.5)

ax[0, 1].legend(loc='upper left', fontsize=7.5, handletextpad=0.1)

# tidy-up and save. 
plt.tight_layout()
plt.savefig('../figures/Fig2_DNA_muts.svg', bbox_inches='tight')
