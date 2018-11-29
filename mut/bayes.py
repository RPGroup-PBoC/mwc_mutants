# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import pickle
import pystan
import bokeh.plotting
from .viz import bokeh_traceplot

class StanModel(object):
    R"""
    Custom StanModel class for crafting and sampling from Stan
    models.
    """
    def __init__(self, file, data_dict=None, samples=None, force_compile=False):
        """
        Parameters
        ----------
        model: str
            Relative path to saved Stan model code. To deter bad habits,
            this class does not accept a string as the model code. Save 
            your Stan models. 
        data_dict: dictionary
            Dictonary of all data block parameters for the model.
        force_compile: bool
            If True, model will be forced to compile. If False, 
            a precompiled file will be loaded if present. 
        """
        if '.pkl' in file:
            s = _load(file)
            self.model = s[0]
            self.samples = s[1]
        else:
            self.model = loadStanModel(file, force=force_compile)
            self.data = data_dict
            self.samples = samples
            self.df = None
        
    def sample(self, data_dict=None, iter=2000, chains=4, return_df=True, **kwargs):
        """
        Samples the assembled model given the supplied data dictionary
        and returns output as a dataframe.
        """
        if data_dict == None:
            data_dict = self.data
        self.chains = chains
        self.iter = iter
        self.samples = self.model.sampling(data_dict, 
                        chains=chains, iter=iter, **kwargs)
        if return_df:
            self.df = self.samples.to_dataframe(diagnostics=True)
            return [self.samples, self.df]
        else:
            return self.samples
    
    # Pickling objects
    def dump(fname):
        """Saves StanFit4Model object and sampling summary as a pickled dictionary."""
        with open(f"{fname.split('.')[0]}.pkl", 'wb') as _file:
            pickle.dump({'model' : self.model, 'fit' : self.samples}, _file, protocol=-1)
                  
    def _load(fname):
        with open(file, 'rb') as _file:
            fit_dict = pickle.load(_file)
        self.model = fit_dict[0]
        self.samples = fit_dict[1]
        return [self.model, self.samples]
    
                  
    def summarize_parameters(self, parnames=[], mass_frac=0.95):
        """
        Summarizes all or a subset of parameters from a Stan model. 
        
        Parameters
        ----------
        parnames: list
            List of desired parnames. If left empty, all parameters 
            are summarized and returned. 
        mass_frac: float [0, 1]
            The probability mass fraction for the HPD. Default is 
            the 95% credible region. 
            
        Returns
        -------
        summary_df: pandas DataFrame
            Dataframe of summarized parameters. The columns are as
            follows:
                parameter = name of parameter in Stan model
                dimension = index (dimension) of the parameter
                mean = mean of samples
                median = median of samples
                mode = parameter value when the log posterior is maximized
                hpd_min = minimum bound of the highest probability density
                    defined by the mass fraction.
                hpd_max = upper bound of the highest probability density
                    defined by the mass fraction
        """
        # Extract the sampling information and find the mode
        samples = self.samples
        fit = samples.extract()
        mode_ind = np.argmax(fit['lp__'])
        
        # Get a list of all parameters defined in the model and assign a dimension
        pars = samples.model_pars
        
        # Convert the dimensions for each parameter to integers. 
        _dims = []
        for d in samples.par_dims:
            if len(d) == 0:
                _dims.append(1)
            else:
                _dims.append(int(d[0]))
    
        par_dims = {p:v for p, v in zip(pars, _dims)}
        if len(parnames) != 0:
            pars = parnames
            desired_pars = {k:v for k, v in par_dims.items() if k in parnames}
            par_dims = desired_pars
        
        # Iterate through each parameter and compute the aggregate properties. 
        df = pd.DataFrame([], columns=['parameter', 'dimension', 'mean'
                                      'mode', 'median', 'hpd_min',
                                      'hpd_max', 'mass_fraction'])          
        for par, dim in par_dims.items():
            par_samples = fit[par]
            if dim == 1:
                par_samples = par_samples[:, np.newaxis]
            for j in range(dim):
                # Compute the summary statistics
                par_mode = par_samples[:, j][mode_ind]
                par_mean = np.mean(par_samples[:, j])
                par_median = np.median(par_samples[:, j])
                hpd_min, hpd_max = compute_hpd(par_samples[:, j], mass_frac=mass_frac)
                
                # Assemble a dictionary to append to the data frame
                par_dict ={'parameter':par,
                          'dimension': j + 1,
                          'mean': par_mean,
                          'mode': par_mode,
                          'median': par_median,
                          'hpd_min': hpd_min,
                          'hpd_max': hpd_max,
                          'mass_fraction': mass_frac}
                df = df.append(par_dict, ignore_index=True)
        df['dimension'] = df['dimension'].astype(int) 
        return df 

        
    # Vizualization    
    def traceplot(self, varnames=None):
        """
        Shows the sampling trace and distributions for desired varnames
        See documentation for mwc.viz.bokeh_traceplot for more details.
        """
        return bokeh_traceplot(self.samples, varnames=varnames) 

     
    def parcoord(self,  varnames=None):
        """
        Creates a parallel coordinate plot with centered parameters. 
        """
        # Get names of unconstrained parameters. 
        parnames = self.samples.unconstrained_param_names()
    
        # Create a new dataframe and rescale the parameters from 0 to 1
        centered_df = pd.DataFrame()
        formatted_parnames = ['lp__']
        for p in parnames:
            par, num = p.split('.')
            formatted_parnames.append(f'{par}[{num}]')
            samp = self.df[f'{par}[{num}]']
            centered_df[f'{par}[{num}]'] = (samp - samp.min()) / (samp.max() - samp.min())
        centered_df['lp__'] = (self.df['lp__'] - self.df['lp__'].min()) / (self.df['lp__'].max() - self.df['lp__'].min())
        centered_df['divergence'] = self.df['divergent__']
        centered_df['samp_idx'] = np.arange(0, len(self.df)) + 1
    

        # Set up the figure canvas. 
        p = bokeh.plotting.figure(width=500, height=200, 
                            x_axis_label='parameter',
                            y_axis_label='normalized',
                            x_range=formatted_parnames)
        
        # Group each sample and plot. 
        for g, d in centered_df.groupby('samp_idx'):
            if d['divergence'].values[0] != 0:
                color ='tomato'
            else:
                color='black'
            p.line(x=formatted_parnames, y=d[formatted_parnames].values[0], line_width=0.1, alpha=0.5, color=color)
        return p

                

def loadStanModel(fname, force=False):
    """Loads a precompiled Stan model. If no compiled model is found, one will be saved."""
    # Identify the model name and directory structure
    rel, sm_dir = fname.split('/stan/')
    sm_name = sm_dir.split('.stan')[0]
    pkl_name = f'{rel}/stan/{sm_name}.pkl' 
    # Check if the model is precompiled
    if (os.path.exists(pkl_name)==True) and (force != True):
        print('Found precompiled model. Loading...')
        model = pickle.load(open(pkl_name, 'rb'))
        print('finished!')
    else:
        print('Precompiled model not found. Compiling model...')
        _path = rel + '/stan/'         
        model = pystan.StanModel(fname, include_paths=_path)        
        print('finished!')
        with open(pkl_name, 'wb') as f:
            pickle.dump(model, f)      
    return model
    

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