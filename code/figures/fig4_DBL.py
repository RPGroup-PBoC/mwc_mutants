# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.viz
import mut.stats
import mut.thermo
mut.viz.plotting_style()
colors = mut.viz.color_selector('mut')
pboc_colors = mut.viz.color_selector('pboc')

# Load the data.
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[(data['class'] == 'WT') & (data['class'] == 'DBL')]
DNA_muts = ['Y20I', 'Q21A', 'Q21M']
IND_muts = ['F164T', 'Q294K', 'Q294V']

# %%
# Set up the figure canvas
fig, ax = plt.subplots(3, 3, figsize=(3.5, 3.5), sharex=True, sharey=True)
axes = ax.ravel()
for a in axes:
    a.set_xlim([1E-8, 1E-2])
    a.set_xscale('log')
    a.set_ylim([-0.15, 1.15])
    a.tick_params(labelsize=8)


for i, dna in enumerate(DNA_muts):
    for j, ind in enumerate(IND_muts):


plt.subplots_adjust(hspace=0.05, wspace=0.05)

