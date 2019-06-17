# -*- coding: utf-8 -*-
# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import mut.bayes
import mut.viz
import mut.thermo
import tqdm
constants = mut.thermo.load_constants()
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()
#%%

# Load the data and  clean
data = pd.read_csv('../../data/csv/RazoMejia2018_data.csv')
data = data[data['repressors']  > 0].copy()
data['repressors'] *= 2
data.rename(columns={'fold_change_A':'fold_change', 'IPTG_uM':'IPTGuM'}, inplace=True)

def infer_empirical_bohr(data, model, groupby=['mutant', 'repressors', 'operator', 'IPTGuM'],
                        verbose=True, force_compile=False, **kwargs):
    """
    Infers the empirical bohr parameter (and relevant correction) for a collection of 
    fold-change measurements
    
    Parameters
    ----------
    data: pandas DataFrame object
        The data from which the empirical bohr will be determined. This should have at least
        a fold-change column and a grouping parameter.
    model: str
        Path to Stan model to load. Model will be compiled if `force_compile`==True.
    groupby: list, optional
        List of identifiers by which to group the supplied data. Default groups by 
        'mutant', 'repressors', 'operator', and 'IPTGuM'
    verbose: bool
        If true, the progress will be printed to screen as a bar. 
    force_compile: bool
        If True, the stan model will be recompiled.
    **kwargs: keyword arguments
        kwargs to be passed to the sampler.
        
    Returns
    -------
    statistics: pandas DataFrame
        Dataframe of statistics for relevant parameters.
    """
    
    # Load the stan model and compile if needed. 
    model = mut.bayes.StanModel(model, force_compile=force_compile)
    
    # Make a storage list for the individual statistics
    fc_stats = []
    
    # Make a quiet or loud iterator. 
    if verbose:
        iter = tqdm.tqdm(data.groupby(groupby))
    else:
        iter = data.groupby(groupby)
        
    # Iter through each grouping and infer
    for g, d in iter: 
        # Define parameters of the reference state
        ref = d['ref_bohr'].unique()[0]
    
        # Assemble the data dictionary and sample the posterior
        data_dict = {'N':len(d),
                     'foldchange': d['fold_change']} 
        fit, samples = model.sample(data_dict, **kwargs)

        # Compute the empirical bohr and delta F 
        samples['empirical_bohr'] = -np.log(samples['fc_mu']**-1 - 1)

        # Identify the extrema
        samples['delta_bohr'] = samples['empirical_bohr'] - ref

        _dbohr_stats = mut.stats.compute_statistics(samples, 
                varnames=['empirical_bohr', 'fc_mu', 'fc_sigma', 'delta_bohr'],  
                logprob_name='lp__')    
        _dbohr_stats['repressors'] = g[0]
        _dbohr_stats['operator'] = g[1]
        _dbohr_stats['IPTGuM'] = g[2]
        # _dbohr_stats['class'] = d['class'].unique()[0]
        fc_stats.append(_dbohr_stats)
    
    return pd.concat(fc_stats)


# %%
# Compute the Delta F relative to R260.
data['ref_bohr'] = mut.thermo.SimpleRepression(R=260, ep_r=data['binding_energy'], 
                    ep_ai=constants['ep_AI'], ka=constants['Ka'], 
                    ki=constants['Ki'], n_sites=constants['n_sites'],
                    effector_conc=data['IPTGuM']).bohr_parameter()

data.head()
R_stats = infer_empirical_bohr(data, '../stan/empirical_F.stan',  groupby=['repressors', 'operator', 'IPTGuM'], 
            **dict(iter=5000, control=dict(adapt_delta=0.99)))
# R_stats.to_csv('../../data/csv/empirical_F_statistics.csv', index=False)

# %%
# Compute the delta F relative to O2
data['ref_bohr'] = mut.thermo.SimpleRepression(R=data['repressors'], ep_r=constants['O2'], 
                    ep_ai=constants['ep_AI'], ka=constants['Ka'], 
                    ki=constants['Ki'], n_sites=constants['n_sites'],
                    effector_conc=data['IPTGuM']).bohr_parameter()

data.head()
O2_stats = infer_empirical_bohr(data, '../stan/empirical_F.stan',  groupby=['repressors', 'operator', 'IPTGuM'], 
            **dict(iter=5000, control=dict(adapt_delta=0.99)))

#%%
# Compute the delta F relative to c = 50 µM 
data['ref_bohr'] = mut.thermo.SimpleRepression(R=data['repressors'], ep_r=data['binding_energy'], 
                    ep_ai=constants['ep_AI'], ka=constants['Ka'], 
                    ki=constants['Ki'], n_sites=constants['n_sites'],
                    effector_conc=50).bohr_parameter()

data.head()
c_stats = infer_empirical_bohr(data, '../stan/empirical_F.stan', groupby=['repressors', 'operator', 'IPTGuM'],
            **dict(iter=5000, control=dict(adapt_delta=0.99)))


#%% Set up the figure canvas
fig, ax = plt.subplots(1, 3, figsize=(7.5, 3))

for a in ax:
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    # a.set_xscale('symlog', linthreshx=1E-2)
    # a.set_xlim([-0.001, 1E4])
    # a.set_xlabel('IPTG [µM]', fontsize=8)
    a.set_ylim([-6, 6])

ax[0].set_ylabel('$\Delta F$ [$k_BT$]', fontsize=8)

# Set the titles to describe the varied parameters
ax[0].set_title('varying $R$ ; $R^{(\mathrm{ref})} = 260$', fontsize=8, 
                backgroundcolor=colors['pale_yellow'])
ax[1].set_title(r'varying $\Delta\varepsilon_{RA}$ ; $\Delta\varepsilon_{RA}^{(\mathrm{ref})} = -13.9\,k_BT$', 
                fontsize=8, backgroundcolor=colors['pale_yellow'])
ax[2].set_title('varying $c$; $c^{(\mathrm{ref})} = 50$ µM', fontsize=8, 
                backgroundcolor=colors['pale_yellow'])

rep_colors = {22: colors['blue'], 60:colors['green'], 124:colors['purple'], 
              260:colors['red'], 1220:colors['dark_brown'], 
              1740:colors['dark_green']}
op_glyphs = {'O1':'s', 'O2':'o', 'O3':'D'}


# ##############################################################################
#  REPRESSOR COPY NUMBER
# ##############################################################################
# Plot the theory curves
rep = np.logspace(0, 4.2)
ax[0].plot(rep / 260, -np.log(rep / 260), color='k', label='theory', lw=1)
ax[0].set_xlim([0.06, 8])
ax[0].set_xscale('log')
ax[0].set_xlabel(r'$\frac{R}{R^{(\mathrm{ref})}}$', fontsize=8)
for g, d in R_stats.groupby(['repressors', 'IPTGuM', 'operator']):
    dF = d[d['parameter']=='delta_bohr']
    fc = d[d['parameter']=='fc_mu']['median'].values[0]
    sig = d[d['parameter']=='fc_sigma']['median'].values[0]
    if (fc > sig) & (1 - fc > sig):
        ax[0].plot(g[0] / 260, dF['median'], marker=op_glyphs[g[-1]], color=rep_colors[g[0]],
        ms=2)
        ax[0].vlines(g[0] / 260, dF['hpd_min'], dF['hpd_max'], lw=1, color=rep_colors[g[0]])


# Plot the operator scaling
dep = np.linspace(-16, -8, 200)
ax[1].plot(dep + 13.9, dep + 13.9, 'k-', lw=1)
for g, d in O2_stats.groupby(['repressors', 'IPTGuM', 'operator']):
    dF = d[d['parameter']=='delta_bohr']
    fc = d[d['parameter']=='fc_mu']['median'].values[0]
    sig = d[d['parameter']=='fc_sigma']['median'].values[0]
    if (fc > sig) & (1 - fc > sig):
        ax[1].plot(constants[g[-1]] + 13.9, dF['median'], marker=op_glyphs[g[-1]], color=rep_colors[g[0]],
        ms=2)
        ax[1].vlines(constants[g[-1]] + 13.9, dF['hpd_min'], dF['hpd_max'], lw=1, color=rep_colors[g[0]])


# Plot the induction curve scaling
c_range = np.logspace(-2, 4, 200)
c_range[0] = 0
pact = mut.thermo.MWC(ka=constants['Ka'], ki=constants['Ki'], ep_ai=constants['ep_AI'],
                      effector_conc=c_range).pact()
pact_ref = mut.thermo.MWC(ka=constants['Ka'], ki=constants['Ki'], ep_ai=constants['ep_AI'],
                      effector_conc=50).pact()

ref_dF = -np.log(pact / pact_ref)
c0 = 50
ax[2].plot(c_range/c0, ref_dF, 'k-', lw=1)
ax[2].set_xscale('log')
ax[2].set_xlim([5E-4, 5E2])

dep = np.linspace(-16, -8, 200)
ax[1].plot(dep + 13.9, dep + 13.9, 'k-', lw=1)
for g, d in c_stats.groupby(['repressors', 'IPTGuM', 'operator']):
    dF = d[d['parameter']=='delta_bohr']
    fc = d[d['parameter']=='fc_mu']['median'].values[0]
    sig = d[d['parameter']=='fc_sigma']['median'].values[0]
    if (fc > sig) & (1 - fc > sig):
        ax[2].plot(g[1]/c0, dF['median'], marker=op_glyphs[g[-1]], color=rep_colors[g[0]],
        ms=2)
        ax[2].vlines(g[1]/c0, dF['hpd_min'], dF['hpd_max'], lw=1, color=rep_colors[g[0]])

#%%
