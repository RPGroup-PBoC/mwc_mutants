# -*- coding: utf-8 -*- import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import mut.thermo
import mut.bayes
import mut.stats
import tqdm
constants = mut.thermo.load_constants()
model = mut.bayes.StanModel('../stan/empirical_deltaF.stan')

# Load the data and isolate to DNA mutants. 
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[data['class'] != 'WT']
data.dropna(inplace=True)

# Define the reference bohr. 
c_0 = 50 # in µM
op_en = [constants[o] for o in data['operator']]
ref = mut.thermo.SimpleRepression(R=data['repressors'],
                                      ep_r=np.array(op_en), ka=constants['Ka'],
                                      ki=constants['Ki'], ep_ai=constants['ep_AI'],
                                      n_sites=constants['n_sites'], effector_conc=c_0)
# Compute the reference strain. 
ref_strain = mut.thermo.SimpleRepression(R=data['repressors'],
                                  ep_r=constants['O2'], ka=constants['Ka'],
                                  ki=constants['Ki'], ep_ai=constants['ep_AI'],
                                  n_sites=constants['n_sites'], 
                                  effector_conc=data['IPTGuM'])

# Add the reference information to the data frame. 
data['F0'] = -ref.bohr_parameter()
data['deltaF_ref'] = data['F0'] + ref_strain.bohr_parameter()

# Set up IPTG idx. 
data['idx'] = data.groupby(['IPTGuM']).ngroup() + 1 
iptg_idx = {i:iptg for i, iptg in zip(data['idx'].unique(), data['IPTGuM'].unique())}
# Set up a storage list for the inference.
df = pd.DataFrame([], columns=['class', 'mutant',
                                'operator', 'repressors',
                                'IPTGuM', 'fc_mu',
                                'fc_mu_min', 'fc_mu_max',
                                'measured_F', 'measured_F_min', 'measured_F_max',
                                'empirical_delF', 'empirical_delF_min',
                                'empirical_delF_max', 'ddF',
                                'ddF_min', 'ddF_max'])


for g, d in tqdm.tqdm(data.groupby(['class', 'mutant', 'repressors', 'operator'])):
    print(g)
    # Set up the data dictionary. 
    data_dict = {'N':len(d),
                'J': len(d['IPTGuM'].unique()),
                'idx':d['idx'], 
                'F0': d['F0'].unique()[0],
                'deltaF_ref': d['deltaF_ref'].unique(),
                'fc': d['fold_change']}
    # Sample the model and convert to data frame. 
    samples, samples_df = model.sample(data_dict=data_dict, iter=10000, thin=10, 
                                      control=dict(adapt_delta=0.95))
    
    for i in range(12):
        # Compute the stats of each:
        fc_median = np.median(samples_df[f'fc_mu[{i+1}]'])
        fc_low, fc_high = mut.stats.compute_hpd(samples_df[f'fc_mu[{i+1}]'], 0.95) 
        F_median = np.median(samples_df[f'measured_F[{i+1}]'])
        F_low, F_high = mut.stats.compute_hpd(samples_df[f'measured_F[{i+1}]'], 0.95)
        empF_median = np.median(samples_df[f'empirical_delF[{i+1}]'])
        empF_low, empF_high = mut.stats.compute_hpd(samples_df[f'empirical_delF[{i+1}]'], 0.95)
        ddf_median = np.median(samples_df[f'deldelF[{i+1}]'])
        ddf_low, ddf_high = mut.stats.compute_hpd(samples_df[f'deldelF[{i+1}]'], 0.95)

        # Assemble the dictionary. 
        samp_dict = {'class':g[0], 'mutant':g[1], 'repressors':g[2], 'IPTGuM':iptg_idx[i+1],
                     'fc_mu':fc_median, 'fc_mu_min':fc_low, 'fc_mu_max':fc_high,
                     'measured_F':F_median, 'measured_F_min':F_low, 'measured_F_max':F_high,
                     'empirical_delF':empF_median, 'empirical_delF_min':empF_low,
                     'empirical_delF_max':empF_high,
                     'ddF': ddf_median, 'ddF_min':ddf_low, 'ddF_max':ddf_high}
        
        df = df.append(samp_dict, ignore_index=True)

