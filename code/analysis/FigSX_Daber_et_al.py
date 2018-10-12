# -*- coding: utf-8 -*-
import sys
import numpy as np
import pystan 
import pandas as pd
sys.path.insert(0, '../../')
import mut.bayes
import mut.stats

# Load the data and assign identifiers
data = pd.read_csv('../../data/csv/Daber2011_data.csv')
DNA_muts = {m: i + 1 for i, m in enumerate(data[data['class']=='DNA']['mutant'].unique())}
IND_muts = {m: i + 1 for i, m in enumerate(data[data['class']=='IND']['mutant'].unique())}
for m, idx in DNA_muts.items():
    data.loc[data['mutant']==m, 'idx'] = idx

for m, idx in IND_muts.items():
    data.loc[data['mutant']==m, 'idx'] = idx

data.loc[data['class']=='WT', 'idx'] = 0
data['idx'] = data['idx'].astype(int)

# Load and compile the stan model. 
model = pystan.StanModel('../stan/Daber_et_al_analysis.stan')

# Separate data into classes. 
dna = data[data['class']=='DNA']
ind = data[data['class']=='IND']
wt = data[data['class']=='WT']

# Assemble the data dictionary
data_dict = dict(J_DNA=len(DNA_muts), J_IND=len(IND_muts), N_DNA=len(dna), N_IND=len(ind),
                N_WT=len(wt), DNA_idx=dna['idx'], IND_idx=ind['idx'], 
                n_sites=2, ep_AI=4.5, c_DNA=dna['IPTGuM'], c_IND=ind['IPTGuM'],
                c_WT=wt['IPTGuM'], DNA_fc=dna['fold_change'], IND_fc=ind['fold_change'],
                WT_fc=wt['fold_change'])

# Sample!
samples = model.sampling(data_dict, iter=50000, thin=100, chains=4)



samples.plot()
