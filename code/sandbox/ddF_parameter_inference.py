# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import tqdm
import mut.thermo
import mut.bayes
constants = mut.thermo.load_constants()

# Load the stan model. 
model = mut.bayes.StanModel('../stan/ddF_parameter_estimation.stan')

# Load the data. 
data = pd.read_csv('../../data/csv/compiled_data.csv')


# Remove the WT. 
data = data[data['class'] != 'WT'].copy()


# Define the reference state. 
c_0 = 50 # in µM
op_en = [constants[op] for op in data['operator']]
data['wt_binding_energy'] = op_en
ref_state = mut.thermo.SimpleRepression(R=data['repressors'], ep_r=data['wt_binding_energy'],
                                       ka=constants['Ka'], ki=constants['Ki'],
                                       n_sites=constants['n_sites'], ep_ai=constants['ep_AI'],
                                       effector_conc=c_0)
ref_strain = mut.thermo.SimpleRepression(R=data['repressors'], ep_r=data['wt_binding_energy'],
                                        ka=constants['Ka'], ki=constants['Ki'],
                                        n_sites=constants['n_sites'], ep_ai=constants['ep_AI'],
                                        effector_conc=data['IPTGuM'])

# compute the empirical ddf 
data['measured_F'] = np.log((1 / data['fold_change']) - 1)
data['empirical_dF'] = -ref_strain.bohr_parameter() - data['measured_F']
data['ddF'] = data['empirical_dF'] + ref_strain.bohr_parameter()

# TODO: Figure out how to deal with this annoying fact. 
# Remove nan values 
data.dropna(inplace=True)

# Instantiate the storage data frame. 
df = pd.DataFrame([], columns=['repressors', 'operator', 'mutant', 'class',
                               'ep_RA_median', 'ep_RA_min', 'ep_RA_max',
                               'Ka_median', 'Ka_min', 'Ka_max',
                               'li_median', 'Ki_min', 'Ki_max',
                               'ep_AI_median', 'ep_AI_min', 'ep_AI_max'])


# Begin the inference. 
for g, d in tqdm.tqdm(data.groupby(['class', 'mutant', 'repressors', 'operator'])):
    print(g) 
    data_dict = {'N':len(d),
            'wt_Ka':constants['Ka'],
            'wt_Ki':constants['Ki'],
            'n_sites':constants['n_sites'],
            'c':d['IPTGuM'],
            'wt_ep_AI': constants['ep_AI'],
            'wt_ep_RA': d['wt_binding_energy'].unique()[0],
            'ddF': d['ddF']}
    
    # Perform sampling and extract. 
    samples, samples_df = model.sample(data_dict)
    
    # Extract parameters 
    parnames = samples.unconstrained_param_names()
    samp_dict = {}
    for i, p in enumerate(parnames):
        samp_dict[f'{p}_median'] = np.median(samples_df[p])
        hpd = mut.stats.compute_hpd(samples_df[p], 0.95)
        samp_dict[f'{p}_min'] = hpd[0]
        samp_dict[f'{p}_max'] = hpd[0]
    samp_dict['class'] = g[0]
    samp_dict['mutant'] = g[1]
    samp_dict['repressors'] = g[2]
    samp_dict['operator'] = g[3]
    df = df.append(samp_dict, ignore_index=True) 
   
df.to_csv('../../data/csv/ddf_full_inference.csv', index=False)
