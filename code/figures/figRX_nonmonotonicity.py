# -*- coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
import bokeh.io
import bokeh.plotting
from bokeh.models import CustomJS
from bokeh.models.widgets import Slider
mut.viz.plotting_style()
bokeh.io.output_notebook()

# %%
# Define the dependent variables
c_range = np.logspace(-2, 4, 500)
kaki_range = np.logspace(-3, 3, 500)

# Define the wild-type parameter values

# Set some initial kas and kis
ka = [10, 100, 500, 1000] # in µM
ki = [0.01, 0.1, 1, 10] # in µM 
epAI = 4.5 # in kT

# Mesh the parameters together. 
c, kaki = np.meshgrid(c_range, kaki_range)

surfs = []
for ka, ki in zip(ka, ki):
    # Define the function for the delta F
    wt_pact = mut.thermo.MWC(ka=ka, ki=ki, ep_ai=epAI, effector_conc=c).pact()
    mut_pact = mut.thermo.MWC(ka=kaki * ka, ki=0.1 * ki, ep_ai=epAI, effector_conc=c).pact()
    delF = np.log(mut_pact / wt_pact)
    surfs.append(delF)
#%% 
# Set up the figure
fig, ax = plt.subplots(3, 4, figsize=(7, 5), sharex=True, sharey=True)
for a in ax.ravel():
    a.xaxis.grid(False)
    a.yaxis.grid(False)
for i in range(4):
    ax[0, i].imshow(surfs[i], vmin=-10, vmax=10, cmap='magma')


#%%
