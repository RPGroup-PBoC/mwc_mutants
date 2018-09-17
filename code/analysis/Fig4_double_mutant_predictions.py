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

# Constants for file naming. 
FIG_NO = 4
OPERATOR = 'O2'
RBS = 'RBS1027'

# %% Load the sampling traces. 
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[data['class']=='DBL']
epRA_samples = pd.read_csv('../../data/csv/Fig2_{}_DNA_binding_energy_samples.csv'.format(OPERATOR))
kaki_samples = pd.read_csv('../../data/csv/Fig3_{}_KaKi_only_samples.csv'.format(OPERATOR))
kaki_epAI_samples = pd.read_csv('../../data/csv/Fig3_O2_KaKi_epAI_samples.csv'.format(OPERATOR))

# Define the sampling constants. 
n_draws = int(1E5)

# Get the mutants
muts = data['mutant'].unique()
IND_muts = []
DNA_muts = []
for m in muts:
    dna, ind = m.split('-')
    IND_muts.append(ind)
    DNA_muts.append(dna)
DNA_muts = list(set(DNA_muts))
IND_muts = list(set(IND_muts))

# Loop through mutants and isolate parameters.
dfs = []
c_range = np.logspace(-2, 4, 200)
for i, dna in enumerate(tqdm.tqdm(DNA_muts)):
    epRA_draws = epRA_samples['ep_RA.{}.{}'.format(dna, constants[RBS])].sample(n_draws, replace=True).values
    for j, ind in enumerate(IND_muts):
        # Determine if ep_AI should be drawn. 
        if ind == 'Q294K':
            _kaki_samps = kaki_epAI_samples.sample(n_draws, replace=True)
            epAI_draws = _kaki_samps['ep_AI.{}'.format(ind)].values
        else: 
            _kaki_samps = kaki_samples.sample(n_draws, replace=True)
            epAI_draws = np.ones(n_draws) * constants['ep_AI']
            
        # Determine Ka/Ki
        ka_draws = _kaki_samps['Ka.{}'.format(ind)].values
        ki_draws = _kaki_samps['Ki.{}'.format(ind)].values
       
        # Compute the fold-change
        cred_region = np.zeros([2, len(c_range)])
        bohr_cred = np.zeros([2, len(c_range)])
        bohr_median = np.zeros(len(c_range))
        for k, c in enumerate(c_range):
            arch = mut.thermo.SimpleRepression(R=constants[RBS], ep_r=epRA_draws,       
                                             ka=ka_draws, ki=ki_draws,ep_ai=epAI_draws,
                                            n_sites=constants['n_sites'], n_ns=constants['Nns'],
                                            effector_conc=c)
            cred_region[:, k] = mut.stats.compute_hpd(arch.fold_change(), 0.95)
            bohr = arch.bohr_parameter()
            bohr_cred[:, k] = mut.stats.compute_hpd(bohr, 0.95)
            bohr_median[k] = np.median(bohr)
            
            # Assemble the dataframe. 
            _df = pd.DataFrame(np.array([c_range, cred_region[0, :], cred_region[1, :],
                                        bohr_cred[0, :], bohr_cred[1, :], bohr_median]).T, 
                                       columns=['IPTGuM', 'hpd_min', 'hpd_max', 
                                                'bohr_min', 'bohr_max', 'bohr_median'])
            _df['DNA_mutant'] = dna
            _df['IND_mutant'] = ind 
            dfs.append(_df)
            
df = pd.concat(dfs, sort=False)
df.to_csv('../../data/csv/Fig{}_{}_double_samples.csv'.format(FIG_NO, OPERATOR))
