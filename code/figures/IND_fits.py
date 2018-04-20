#%%
# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.thermo
import mut.stats
import mut.viz
colors = mut.viz.pub_style()

# Load the data and restrict to inducer mutants.
data = pd.read_csv('../../data/csv/compiled_data.csv')
IND = data[(data['class']=='IND') | (data['class'] == 'WT')]

# Load the fit statistics and the flat chains.
samples = pd.read_csv('../../data/mcmc/IND_O2_KaKi_fit_chains.csv')
stats = pd.read_csv('../../data/mcmc/IND_O2_KaKi_fit_statistics.csv')