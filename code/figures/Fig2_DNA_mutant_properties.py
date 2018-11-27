# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.thermo
import mut.viz
constants = mut.thermo.load_constants()
colors = mut.viz.color_selector('mut')
pboc = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load data and sampling statistics
data = pd.read_csv('../../data/csv/summarized_data.csv')
data = data[((data['class']=='WT') | (data['class']=='DNA')) & (data['operator']=='O2')]
stats = pd.read_csv('../../data/csv/Fig2_O2_DNA_binding_energy_stats.csv')

# Define IPTG concentration range.
c_range = np.logspace(-2, 4, 200)
r_kdna_range = np.logspace(-2, 6, 200)

# Define functions for calculation of collapse curves. 
def leakiness(R_Kdna, ep_ai=constants['ep_AI']):
    return (1 + (1 / (1 + np.exp(-ep_ai))) * R_Kdna)**-1 
    
def saturation(R_Kdna, Ka=constants['Ka'], Ki=constants['Ki'],
               ep_ai=constants['ep_AI'], n_sites=constants['n_sites']):
    return (1 + (1 / (1 + np.exp(-ep_ai)*(Ka/Ki)**n_sites)) * R_Kdna)**-1
    
def dynamic_range(R_Kdna, ep_ai=constants['ep_AI'], Ka=constants['Ka'],
                 Ki=constants['Ki'], n_sites=constants['n_sites']):
    return saturation(R_Kdna, Ka, Ki, ep_ai, n_sites) - leakiness(R_Kdna, ep_ai)

def ec50(R_Kdna, Ka=constants['Ka'], Ki=constants['Ki'], n_sites=constants['n_sites'],
        ep_ai=constants['ep_AI']):
    
    # Break it into pieces
    repression = 1 + R_Kdna
    numer = repression + (Ka/Ki)**n_sites * (2 * np.exp(-ep_ai) + repression)
    denom = 2 * repression + np.exp(-ep_ai) + (Ka/Ki)**n_sites * np.exp(-ep_ai)
    
    # Assemble
    ec50_numer =  (Ka/Ki) - 1
    ec50_denom = (Ka/Ki) - (numer / denom)**(1 / n_sites)


    return Ka * ((ec50_numer / ec50_denom) - 1)

def effective_hill(R_Kdna, Ka=constants['Ka'], Ki=constants['Ki'], ep_ai=constants['ep_AI'],
                  n_sites=constants['n_sites']):
    # Break it into pieces
    c = ec50(R_Kdna, Ka=Ka, Ki=Ki, n_sites=n_sites, ep_ai=ep_ai)
    pact = (1 + (c / Ka))**n_sites / ((1 + c/Ka)**n_sites + np.exp(-ep_ai) * (1 + c / Ki)**n_sites) 
    fc = (1 + pact * R_Kdna)**-1
    leak = leakiness(R_Kdna, ep_ai=ep_ai)
    expanded_ka = (1 + c / Ka)
    expanded_ki = (1 + c / Ki)
    
    # Assemble the parts
    prefactor = -fc**2 * R_Kdna * 2 * c * np.exp(-ep_ai)
    numer = (1 / Ka) * expanded_ka * expanded_ki**2 - (1 / Ki) * expanded_ka**2 * expanded_ki
    denom = (expanded_ka**2 + np.exp(-ep_ai) * expanded_ki**2)**2
    return (2 / (fc - leak)) * prefactor * (numer / denom)
    
    
# Instantiate the figure. 
fig, ax = plt.subplots(3, 3, figsize=(7.5, 7.5))

# Format axes
titles = ['(A)', '(B)', '(C)', '(D)', '(E)', '(F)', '(G)', '(H)', '(I)']
for i, a in enumerate(ax.ravel()):
    a.xaxis.set_tick_params(labelsize=6.5)
    a.yaxis.set_tick_params(labelsize=6.5)
    a.set_xscale('log')
    a.text(-0.45, 1.02, f'{titles[i]}', fontsize=6.5, transform=a.transAxes)
    if i > 2:
        a.set_xlim([1E-2,1E4])
        a.set_xlabel(r'$\frac{R}{N_{NS}}e^{-\beta\Delta\varepsilon_{RA}}$', fontsize=8)
    else:
        a.set_xlim([1E-8, 1E-2])
        a.set_xlabel('IPTG [M]', fontsize=6.5)
        a.set_ylabel('fold-change', fontsize=6.5)
ax[0, 0].axis('off')
ax[0, 2].axis('off')
# Conditional formatting
ax[0, 1].set_xscale('log')
ax[1, 0].set_xscale('log')
ax[2, 0].set_yscale('log')

#Add appropriate labels. 
ax[1, 0].set_ylabel('leakiness', fontsize=6.5)
ax[1, 1].set_ylabel('saturation', fontsize=6.5)
ax[1, 2].set_ylabel('dynamic range', fontsize=6.5)
ax[2, 0].set_ylabel('[EC$_{50}$] [µM]', fontsize=6.5)
ax[2, 1].set_ylabel('effective Hill coefficient', fontsize=6.5)

# #########################################
# INDUCTION PROFILES
# #########################################
c_range = np.logspace(-8, -2, 200)

# Plot the O2 R260 induction profiles. 
for g, d in data[data['repressors']==260].groupby(['mutant']):
    ax[0, 1].errorbar(d['IPTGuM'] / 1E6, d['mean'], d['sem'], 
                      fmt='o', ms=2, markerfacecolor=colors[g.upper()],
                     linestyle='none', label=g.upper(),
                     color=colors[g.upper()], lw=0.5)
    
    # Compute and plot the fit line
    if g != 'wt':
        epRA = np.array(list(stats[stats['parameter']==f'ep_RA.{g.upper()}.260'].values[0][1:]))
        
    else:
        epRA = np.array([-13.9, -13.9, -13.9])
        
    epRA_mesh, c_mesh = np.meshgrid(epRA, c_range)
    arch = mut.thermo.SimpleRepression(R=260, ep_r=epRA_mesh, effector_conc=c_mesh,
                                      ep_ai=constants['ep_AI'], n_ns=constants['Nns'],
                                      ka=constants['Ka']/1e6, ki=constants['Ki']/1e6,
                                      n_sites=constants['n_sites']).fold_change()
    ax[0, 1].plot(c_range, arch[:, 0 ], lw=1, color=colors[g.upper()], label='__nolegend__')
    ax[0, 1].fill_between(c_range, arch[:, 1], arch[:, 2], alpha=0.5, color=colors[g.upper()], label='__nolegend__')

# Define the property-axis mapping
ax_map = {'leakiness':ax[1, 0], 'saturation':ax[1, 1],
         'dynamic_range':ax[1,2], 'EC50':ax[2, 0], 'effective_hill':ax[2, 1]}
_glyphs = ['o', 's', 'D', '^']
glyphs = {r:g for r, g in zip(data['repressors'].unique(), _glyphs)}


# Add a legend on plot 1 regarding glyph shapes
for r, g in glyphs.items():
    ax[1,0].plot([], [], g, label=int(r), ms=3, color='slategray')
ax[0, 1].legend(fontsize=6.5, handletextpad=0.05, loc='upper left', ncol=1,
               columnspacing=0.02)
_leg = ax[1, 0].legend(fontsize=6.5, loc='upper right', handletextpad=0.1, 
                      title='rep. / cell')
_leg.get_title().set_fontsize(8)

# #########################################
# PLOT AND COMPUTE PROPERTIES
# #########################################
# Compute and plot the properties.
_ = ax[1, 0].plot(r_kdna_range, leakiness(r_kdna_range), 'k-', lw=1)
_ = ax[1, 1].plot(r_kdna_range, saturation(r_kdna_range), 'k-', lw=1)
_ = ax[1, 2].plot(r_kdna_range, dynamic_range(r_kdna_range), 'k-', lw=1)
_ = ax[2, 0].plot(r_kdna_range, ec50(r_kdna_range), 'k-', lw=1)
_ = ax[2, 1].plot(r_kdna_range, effective_hill(r_kdna_range), 'k-', lw=1)

# Go through each mutant and compute the properties. 
for g, d in data.groupby(['mutant', 'repressors']):
    if g[0] == 'wt':
        ep_RA = np.array([-13.9, -13.9, -13.9])
        ep_RA_fit = np.array([-13.9, -13.9, -13.9])
    else:
        ep_RA = stats[stats['parameter']==f'ep_RA.{g[0]}.{260}'][['mode', 'hpd_min', 'hpd_max']].values[0]
        ep_RA_fit = stats[stats['parameter']==f'ep_RA.{g[0]}.{int(g[1])}'][['mode', 'hpd_min', 'hpd_max']].values[0]
         
    # Compute R/K_dna
    R_Kdna = (g[1] / constants['Nns']) * np.exp(-ep_RA)
     
    # Plot the properties. 
    if g[1] == 260:
        label=g[0].upper()
        face = 'w'
    else:
        label='__nolegend__'
        face = colors[g[0]]
   
    # Leakiness
    ax[1, 0].errorbar(R_Kdna[0], d[d['IPTGuM']==0]['mean'], d[d['IPTGuM']==0]['sem'], lw=1, fmt=glyphs[g[1]], color=colors[g[0].upper()],
                     markerfacecolor=face, markeredgecolor=colors[g[0].upper()], linestyle='none', markersize=3,
                     label=label)
    ax[1, 0].hlines(d[d['IPTGuM']==0]['mean'], R_Kdna[2], R_Kdna[1], lw=1, color=colors[g[0].upper()], label='__nolegend__')
    
    # Saturation
    ax[1, 1].errorbar(R_Kdna[0], d[d['IPTGuM']==5000]['mean'], d[d['IPTGuM']==5000]['sem'], lw=1, fmt=glyphs[g[1]], color=colors[g[0].upper()],
                     markerfacecolor=face, markeredgecolor=colors[g[0].upper()], linestyle='none', markersize=3) 
    ax[1, 1].hlines(d[d['IPTGuM']==5000]['mean'], R_Kdna[1], R_Kdna[2], lw=1, color=colors[g[0].upper()])
   
    # Dynamic Range
    _ = ax[1, 2].plot(R_Kdna[0], d[d['IPTGuM']==5000]['mean'].values - d[d['IPTGuM']==0]['mean'].values,
                     marker=glyphs[g[1]], markersize=3, markerfacecolor=face, markeredgecolor=colors[g[0].upper()])

    # Plot the inferred parameters.
    arch = mut.thermo.SimpleRepression(R=g[1], ep_r=ep_RA_fit, ka=constants['Ka'], ki=constants['Ki'],
                                      n_sites=constants['n_sites'], n_ns=constants['Nns'], 
                                       ep_ai=constants['ep_AI'], effector_conc=0).compute_properties()
    
    # EC50
    _ = ax[2, 0].plot(R_Kdna[0], arch['EC50'][0], marker=glyphs[g[1]], markerfacecolor=face, markeredgecolor=colors[g[0].upper()],
                     markersize=3)
    _ = ax[2, 0].hlines(arch['EC50'][0], R_Kdna[1], R_Kdna[2], lw=1, color=colors[g[0].upper()])
    _ = ax[2, 0].vlines(R_Kdna[0], arch['EC50'][1], arch['EC50'][2], lw=1, color=colors[g[0].upper()])
    
    
    # Effective Hill
    _ = ax[2, 1].plot(R_Kdna[0], arch['effective_hill'][0], marker=glyphs[g[1]], markerfacecolor=face, markeredgecolor=colors[g[0].upper()],
                     markersize=3)
    _ = ax[2, 1].hlines(arch['effective_hill'][0], R_Kdna[1], R_Kdna[2], lw=1, color=colors[g[0].upper()])
    _ = ax[2, 1].vlines(R_Kdna[0], arch['effective_hill'][1], arch['effective_hill'][2], lw=1, color=colors[g[0].upper()])
     
# #########################################
# DATA COLLAPSE
# #########################################
ax[2, 2].set_xscale('linear')
ax[2, 2].set_xlim([-8, 8])
ax[2, 2].set_xlabel('Bohr parameter [$k_BT$]', fontsize=6.5)
ax[2, 2].set_ylabel('fold-change', fontsize=6.5)

# Iterate through each mutant and compute the bohr parameter. 
for g, d in data.groupby(['mutant', 'repressors']):
    # Isoalte the DNA binding energy
    if g[0] != 'wt':
        epRA = np.array(list(stats[stats['parameter']==f'ep_RA.{g[0]}.260'].values[0][1:]))
    else:
        epRA = np.array([-13.9, -13.9, -13.9])
    epRA_mesh, c_mesh = np.meshgrid(epRA, d['IPTGuM'])
    
    bohr = mut.thermo.SimpleRepression(R=g[1], ep_r=epRA_mesh, effector_conc=c_mesh,
                                      ep_ai=constants['ep_AI'], n_sites=constants['n_sites'],
                                      n_ns=constants['Nns'], ka=constants['Ka'], ki=constants['Ki']).bohr_parameter()
    if g[1] == 260:
        face='w'
    else:
        face=colors[g[0].upper()]

    # Plot the bohr parameter. 
    _ = ax[2, 2].errorbar(bohr[:, 0], d['mean'], d['sem'], color=colors[g[0].upper()], linestyle='none',
                        markerfacecolor=face, markeredgecolor=colors[g[0].upper()], fmt=glyphs[g[1]], 
                        label='__nolegend__', ms=3, markeredgewidth=1)
    _ = ax[2, 2].hlines(d['mean'], bohr[:, 1], bohr[:, 2], lw=1, color=colors[g[0].upper()])
    
 # Plot the master curve for collapse
bohr_range = np.linspace(-8, 8, 200)
collapse = (1 + np.exp(-bohr_range))**-1                                                                  
ax[2, 2].plot(bohr_range, collapse, 'k-', lw=1, label='__nolegend__', zorder=2000)
# #########################################
# SAVING AND FORMATTING
# #########################################
plt.subplots_adjust(hspace=0.5, wspace=0.6)
plt.savefig('Fig1_properties_collapse.svg')


