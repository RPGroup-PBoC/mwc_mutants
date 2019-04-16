# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import mut.viz
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load the posterior check data. 
data = pd.read_csv('../../data/csv/IND_posterior_predictive_samples.csv')


# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig = plt.figure(figsize=(7, 3))
gs = gridspec.GridSpec(6, 10)
ax1 = fig.add_subplot(gs[0:2, 0:2])
ax2 = fig.add_subplot(gs[2:4, 0:2])
ax3 = fig.add_subplot(gs[2:4, 2:4])
ax4 = fig.add_subplot(gs[4:6, 0:2])
ax5 = fig.add_subplot(gs[4:6, 2:4])
ax6 = fig.add_subplot(gs[4:6, 4:6])
ax7 = fig.add_subplot(gs [:, 6:])
ax = [ax1, ax2, ax3, ax4, ax5, ax6, ax7]

# Apply special formatting where needed
for a in ax:
    a.xaxis.set_tick_params(labelsize=6)
    a.yaxis.set_tick_params(labelsize=6)

for a in [ax1, ax3, ax5, ax6]:
    a.set_yticks([])
    a.set_xticklabels([])
ax2.set_xlabel('$K_A$ [µM]', fontsize=8)
ax2.set_ylabel('$K_I$ [µM]', fontsize=8)
ax3.set_xlabel('$K_I$ [µM]', fontsize=8)
plt.tight_layout()


