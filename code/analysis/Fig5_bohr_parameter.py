# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import pystan
sys.path.insert(0, '../../')
import mut.bayes
import mut.stats
import tqdm
constants = mut.thermo.load_constants()

# Define the figure constants
FIG_NO = 5
OPERATOR = 'O2'
RBS = 'RBS1027'
n_draws = int(1E6)

# Load the raw data
data = pd.read_csv('../../data/csv/compiled_data.csv')
muts = data[data['class']=='DBL']['mutant'].unique()

# Load the samples
epRA_samples = pd.read_csv('../../data/csv/Fig2_O2_DNA_binding_energy_samples.csv')
epRA_stats = pd.read_csv('../../data/csv/Fig2_O2_DNA_binding_energy_stats.csv')
kaki_samples = pd.read_csv('../../data/csv/Fig3_O2_KaKi_only_samples.csv')
kaki_stats = pd.read_csv('../../data/csv/Fig3_O2_KaKi_only_stats.csv')
kaki_epAI_samples = pd.read_csv('../../data/csv/Fig3_O2_KaKi_epAI_samples.csv')
kaki_epAI_stats = pd.read_csv('../../data/csv/Fig3_O2_KaKi_epAI_stats.csv')

# Identify the double mutants. 
dfs = []
DNA_muts = []
IND_muts = []
for m in muts:
    dna, ind = m.split('-')
    DNA_muts.append(dna)
    IND_muts.append(ind)
DNA_muts = list(set(DNA_muts))
IND_muts = list(set(IND_muts))

bohr_df = pd.DataFrame([])
for i, dna in enumerate(DNA_muts):
    epRA_samps = epRA_samples['ep_RA.{}.{}'.format(dna, constants[RBS])].sample(n_draws, replace=True).values
    for j, ind in enumerate(IND_muts):
        if ind == 'Q294K':
            _kaki_samps = kaki_epAI_samples.sample(n_draws, replace=True)
            epAI_samps = _kaki_samps['ep_AI.{}'.format(ind)].values
            Ka_mode = kaki_epAI_stats[kaki_epAI_stats['parameter']=='Ka.{}'.format(ind)]['mode'].values[0]
            Ki_mode = kaki_epAI_stats[kaki_epAI_stats['parameter']=='Ki.{}'.format(ind)]['mode'].values[0]
            epAI_mode = kaki_epAI_stats[kaki_epAI_stats['parameter']=='ep_AI.{}'.format(ind)]['mode'].values[0]
            
        else:
            Ka_mode = kaki_stats[kaki_stats['parameter']=='Ka.{}'.format(ind)]['mode'].values[0]
            Ki_mode = kaki_stats[kaki_stats['parameter']=='Ki.{}'.format(ind)]['mode'].values[0] 
            epAI_mode = constants['ep_AI']
            _kaki_samps = kaki_samples.sample(n_draws, replace=True)
            epAI_samps = np.ones(n_draws) * constants['ep_AI']
        ka_draws = _kaki_samps['Ka.{}'.format(ind)].values
        ki_draws = _kaki_samps['Ki.{}'.format(ind)].values
        
        # Compute the most-likely bohr parameter value. 
        
        bohr_mode = mut.thermo.SimpleRepression(R=constants['RBS1027'], 
                            ep_r=epRA_stats[epRA_stats['parameter']=='ep_RA.{}.260'.format(dna)]['mode'].values[0],
                            ka=Ka_mode, ki=Ki_mode, ep_ai=epAI_mode, n_sites=constants['n_sites'],
                            n_ns=constants['Nns'], effector_conc=data['IPTGuM'].unique()).bohr_parameter()
        # Compute the delta bohr for each draw
        for k, c in enumerate(data['IPTGuM'].unique()):
            bohr = mut.thermo.SimpleRepression(R=constants['RBS1027'], ep_r=epRA_samps,
                                              effector_conc=c, ep_ai=epAI_samps,
                                              ka=ka_draws, ki=ki_draws,
                                              n_sites=constants['n_sites'], 
                                               n_ns=constants['Nns']).bohr_parameter()
            cred_region = mut.stats.compute_hpd(bohr, mass_frac=0.95)

            bohr_df = bohr_df.append({'mode': bohr_mode[k], 'hpd_min':cred_region[0], 'hpd_max':cred_region[1], 'mutant':'{}-{}'.format(dna, ind), 'IPTGuM':c},
                                ignore_index=True)
bohr_df.to_csv('../../data/csv/Fig{}_{}_DBL_predicted_Bohr.csv'.format(FIG_NO, OPERATOR)) 

        
        
        
        