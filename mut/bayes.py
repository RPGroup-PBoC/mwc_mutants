# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import scipy.special
import scipy.optimize
import statsmodels.tools.numdiff as smnd

def chains_to_dataframe(fit, varnames=None):
    """
    Converts the generated traces from MCMC sampling to a tidy
    pandas DataFrame.

    Parameters
    ----------
    fit : pystan sampling output
        The raw MCMC output.
    var_names : list of str
        Names of desired parameters. If `None`, all parameters will be
        returned.

    Returns
    -------
    df : pandas DataFrame
        Pandas DataFrame containing all samples from the MCMC.
    """

    data = fit.extract()
    keys = list(data.keys())
    if varnames == None:
        varnames = [k for k in keys if 'lp__' not in k]
    else:
        varnames = fit.unconstrained_param_names()

    samples = {}
    for i, key in enumerate(varnames):
        # Get the shape.
        dim = np.shape(data[key])
        if len(dim) == 2:
            for j in range(dim[-1]):
                samples['{}.{}'.format(key, j+1)] = data[key][:, j]
    
        else:
            samples[key] = data[key]
    samples['logp'] = data['lp__']
    return pd.DataFrame(samples)
    
def assemble_StanModelCode(model_file, function_file):
    """
    Returns a string of the stan model code from a model and function file

    Parameters 
    -----------
    model_file: str
        Path to the model file.
    function_file: str
        Path to the file of functions.

    Returns
    -------
    model_code: str
        String of the stitched together model code.
    """
    lines = []
    files = [function_file, model_file]
    for f in files:
        with open(f, 'r') as file:
            out = file.read().splitlines()
            for line in out:
                lines.append(line) 
    model_code = """\n"""
    for line in lines:
        model_code += line + '\n'
    return model_code


def longform_mcmc_df(dataframe, split_pattern='.', root='logp',
                    idx_name='mutant', param_name='parameter'):
    """
    Convert a dataframe of MCMC samples to longform tidy format.
    
    Parameters
    ----------
    dataframe: pandas DataFrame object.
        The data frame of MCMC samples 
    
    """
    
    # Perform initial melting
    melt = dataframe.melt([root])
    
    # Split the variables given identifier
    idx = [v.split(split_pattern)[1] for v in melt['variable']]
    param = [v.split(split_pattern)[0] for v in melt['variable']]
    
    # Assign the splits
    melt[idx_name] = idx
    melt[param_name] = param
    
    # Drop unnecessary column and perform final melting operation
    melt.drop('variable', axis=1, inplace=True)
    longform = melt.melt([root, idx_name, param_name])
    return longform 