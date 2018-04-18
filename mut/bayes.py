# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import pymc3 as pm
import scipy.special
import scipy.optimize
import statsmodels.tools.numdiff as smnd
import theano.tensor as tt

def chains_to_dataframe(fit, var_names=None):
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
    varnames = [k for k in keys if 'lp__' not in k]
    samples = {}
    for i, key in enumerate(varnames):
        # Get the shape.
        dim = np.shape(data[key])
        if len(dim) == 2:
            for j in range(dim[-1]):
                samples['{}.{}'.format(key, j+1)] = data[key][:, j]
    
        else:
            samples[key] = data[key]
            
    # compute the log_post. 
    new_keys = samples.keys()
    logp = []
    for j in range(dim[0]):
        _samples = [samples[k][j] for k in fit.unconstrained_param_names()]
        if (np.nan not in _samples) & (np.inf not in samples):
            logp.append(fit.log_prob([samples[k][j] for k in fit.unconstrained_param_names()]))

    samples['logp'] = logp
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

