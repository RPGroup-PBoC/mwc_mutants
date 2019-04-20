# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np 
import pandas as pd 
import mut.stats 
import mut.bayes
import mut.thermo
import tqdm
constants = mut.thermo.load_constants()

# Load the compiled data
data = pd.read_csv('../../data/csv/compiled_data.csv')
data.dropna(inplace=True)

# Compute the reference bohr. 
ops = [constants[op] for op in data['operator']]
wt_bohr = mut.thermo.SimpleRepression(R=data['repressors'], ep_r=ops, 
                                       ka=constants['Ka'], ki=constants['Ki'],
                                       ep_ai=constants['ep_AI'], 
                                       effector_conc=data['IPTGuM']).bohr_parameter()
data['ref_bohr'] = wt_bohr

# Assign unique identifiers. 
idx = data.groupby(['mutant', 'repressors', 'operator', 'IPTGuM']).ngroup() + 1 
data['idx'] = idx
data.sort_values('idx', inplace=True)
samples_dfs = []
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
        fc_ref = (1 + np.exp(-ref))**-1
    
        # Assemble the data dictionary and sample the posterior
        data_dict = {'N':len(d),
                     'foldchange': d['fold_change']} 
        fit, samples = model.sample(data_dict, **kwargs)


        # Sample fake data around the mean with the average sigma to estimate error
        sim_data = np.random.normal(fc_ref, samples['fc_sigma'].median(), 
                              len(d))
        _, samples_sim = model.sample(dict(N=len(sim_data),foldchange=sim_data),
                         control=dict(adapt_delta=0.95))

        # Compute the empirical bohr and delta F 
        samples['sim_empirical_bohr'] = -np.log(samples_sim['fc_mu']**-1 -1)
        samples['empirical_bohr'] = -np.log(samples['fc_mu']**-1 - 1)
        samples['correction'] = ref - samples['sim_empirical_bohr'].median()
        # Identify the extrema
        # extrema = (samples['fc_mu'] < samples['fc_sigma']).astype(int) + (1-samples['fc_mu'] < samples['fc_sigma']).astype(int)

        # # Compute the empirical bohr parameter and the delta bohr
        # samples['empirical_bohr'] = -np.log((samples['fc_mu'])**-1 - 1)
        # samples['delta_bohr'] = ref - samples['empirical_bohr']

        # # Compute the delta F error of the reference, given the sigma
        # delta_F_ref_upper = np.nan_to_num(ref + np.log((fc_ref + samples['fc_sigma'])**-1 - 1))
        # delta_F_ref_lower = np.nan_to_num(ref + np.log((fc_ref -  samples['fc_sigma'])**-1 - 1))
        # samples['correction'] = (delta_F_ref_upper-delta_F_ref_lower) * extrema
        samples['delta_bohr_corrected'] = ref - samples['empirical_bohr'] 
        samples['delta_bohr_corrected2'] = samples['sim_empirical_bohr'] - samples['empirical_bohr']

        _dbohr_stats = mut.stats.compute_statistics(samples, varnames=['empirical_bohr', 'fc_mu', 
                                                                       'fc_sigma', 'delta_bohr_corrected',
                                                                       'delta_bohr_corrected2', 'correction'], 
                                               logprob_name='lp__')    
        _dbohr_stats['mutant'] = g[0]
        _dbohr_stats['repressors'] = g[1]
        _dbohr_stats['operator'] = g[2]
        _dbohr_stats['IPTGuM'] = g[3]
        _dbohr_stats['class'] = d['class'].unique()[0]
        fc_stats.append(_dbohr_stats)
    
    return pd.concat(fc_stats)


fc_stats = infer_empirical_bohr(data, '../stan/empirical_F.stan',  **dict(iter=1000, control=dict(adapt_delta=0.99)))
fc_stats.to_csv('../../data/csv/empirical_F_statistics.csv', index=False)

