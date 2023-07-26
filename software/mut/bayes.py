# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import pickle
import pystan
import tqdm
from .stats import compute_statistics

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
