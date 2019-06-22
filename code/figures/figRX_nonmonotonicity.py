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
colors = mut.viz.color_selector('pboc')

# %%
# Define the dependent variables
c_range = np.logspace(-2, 4, 500)
kaki_range = np.logspace(-3, 3, 500)

# Define the wild-type parameter values

# Compute the components of the delta F as a function of c
ka = 200
ki = 1
ep_ai = 4.5
theta1 = 0.1
theta2 = 4
theta3 = 4
theta4 = 2
theta5 = 0.01
theta6=0.1

def dF_dc(ka, ki, ep_ai, c):
    numer = (ka - ki) * (c + ki)
    denom = (c + ka) * ((c + ka)**2 * ki**2 + np.exp(-ep_ai) * ka**2 * (c + ki)**2)
    return numer / denom

def F(ka, ki, ep_ai, c):
    numer = (1 + c/ka)**2
    denom = numer + np.exp(-ep_ai) * (1 + c/ki)**2
    return - np.log(numer/denom)

# %% 
fig, ax = plt.subplots(2, 3, sharex=True, figsize=(7.5, 4.5))
for a in ax.ravel():
    a.set_xscale('log')
    a.set_xlabel('$c \, / \, K_A^\mathrm{(wt)}$', fontsize=12)

wt = 2 * np.exp(-ep_ai) * ka**2 * dF_dc(ka ,ki, ep_ai, c_range)
mut1 = 2 * np.exp(-ep_ai) * ka**2 * dF_dc(ka*theta1 ,ki*theta1, ep_ai, c_range)
mut2 = 2 * np.exp(-ep_ai) * ka**2 * dF_dc(ka*theta2 ,ki*theta2, ep_ai, c_range)
mut3 = 2 * np.exp(-ep_ai) * ka**2 * dF_dc(ka*theta3 ,ki*theta4, ep_ai, c_range)
mut4 = 2 * np.exp(-ep_ai) * ka**2 * dF_dc(ka*theta5 ,ki*theta6, ep_ai, c_range)

Fwt = F(ka, ki, ep_ai, c_range)
Fmut1 = F(ka*theta1, ki*theta1, ep_ai, c_range)
Fmut2 = F(ka*theta2, ki*theta2, ep_ai, c_range)
Fmut3 = F(ka*theta3, ki*theta4, ep_ai, c_range)
Fmut4 = F(ka*theta5, ki*theta6, ep_ai, c_range)

ax[0, 1].plot(c_range / ka, Fwt, label='wild type')
ax[0, 1].plot(c_range / ka, Fmut1, label=r'$\theta=0.1$')
ax[0, 1].plot(c_range / ka, Fmut2, label=r'$\theta=4$')
ax[0, 2].plot(c_range / ka, Fwt - Fwt)
ax[0, 2].plot(c_range / ka, Fmut1 - Fwt)
ax[0, 2].plot(c_range / ka, Fmut2 - Fwt)

ax[1, 1].plot(c_range / ka, Fwt, label='wild type')
ax[1, 1].plot(c_range / ka, Fmut3, label=r'2 $K_A^\mathrm{(wt)} / K_I^\mathrm{(wt)}$')
ax[1, 1].plot(c_range / ka, Fmut4, label=r'0.1 $K_A^\mathrm{(wt)} / K_I^\mathrm{(wt)}$')

ax[1, 2].plot(c_range / ka, Fwt - Fwt)
ax[1, 2].plot(c_range / ka, Fmut3 - Fwt)
ax[1, 2].plot(c_range / ka, Fmut4 - Fwt)

for i in range(2):
    ax[i, 0].axis(False)
    ax[i, 1].set_ylabel('$F\, / \, k_BT$')
    ax[i, 2].set_ylabel('$\Delta F \, / \, k_BT$')

ax[1, 0].set_title('different $K_A/K_I$', backgroundcolor=colors['pale_yellow'], y=0.8)
ax[0, 0].set_title('same $K_A/K_I$', backgroundcolor=colors['pale_yellow'], y=0.8)
ax[0, 1].legend(bbox_to_anchor=(-0.55, 0.8), handlelength=0.5)
ax[1, 1].legend(bbox_to_anchor=(-0.45, 0.8), handlelength=0.5)

# Add panel labels
fig.text(0.1, 0.85, '(A)', fontsize=10)
fig.text(0.1, 0.42, '(B)', fontsize=10)
plt.subplots_adjust(wspace=0.4)
plt.savefig('../../figures/FigR1_kaki_monotonicity.pdf', bbbox_inches='tight')

#%%


#%%
