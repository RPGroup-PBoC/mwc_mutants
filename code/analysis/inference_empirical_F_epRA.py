# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import mut.bayes
import mut.stats
import mut.thermo
constants = mut.thermo.load_constants()

# Load the data
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[data['class']=='DNA'].copy()
data.dropna(inplace=True)

# Load the stan model. 
model = mut.bayes.StanModel('../stan/empirical_F_epRA_estimation.stan')

# Isolate to the correct fitting strain
fit_strain = data[(data['repressors']==260) & 
                  (data['operator']=='O2')].copy()

# Compute the reference state and insert into the data frame. 
ep_r = np.array([constants[op] for op in fit_strain['operator'].values])
arch = mut.thermo.SimpleRepression(R=fit_strain['repressors'], ep_r=ep_r,
                                  ka=constants['Ka'], ki=constants['Ki'],
                                  n_sites=constants['n_sites'], 
                                  ep_ai=constants['ep_AI'],
                                  effector_conc=fit_strain['IPTGuM'])
fit_strain['ref_bohr'] = arch.bohr_parameter()
fit_strain['ref_foldchange'] = arch.fold_change()

# Add identifiers to concentrations, independent of mutation. 
fit_strain['idx'] = fit_strain.groupby(['IPTGuM']).ngroup() + 1

# Assemble the data dictionary. 
_x = fit_strain[fit_strain['mutant']=='Q21M']
data_dict = {'N':len(_x), 'J':_x['idx'].max(),
            'idx':_x['idx'], 'ref_bohr':_x['ref_bohr'],
            'ref_foldchange':_x['ref_foldchange'],
            'foldchange':_x['fold_change']}
fit, samples = model.sample(data_dict)
