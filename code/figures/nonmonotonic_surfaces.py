# -*- coding: utf-8 -*-
# %%
import numpy as np
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
mut.viz.plotting_style()
colors = mut.viz.color_selector('pboc')

ka = 200 
ki = 1 
ep_ai = 5

# Define functions for wild-type and initial cases
def pact(c, ka, ki, ep_ai):
    numer = (1 + c/ka)**2
    denom = numer + np.exp(-ep_ai) * (1 + c/ki)**2
    return numer/denom

def delF(c, theta1, theta2):
    pact_mut = pact(c, theta1 * ka, theta2 * ki, ep_ai)
    pact_wt = pact(c, 139, 0.53, 4.5)
    return - np.log(pact_mut / pact_wt)

def derivF(c, theta1, theta2):
    ka = 139
    ka_mut = theta1 * ka
    ki = 0.53
    ki_mut = theta2 * ki
    ep_ai = 4.5
    wt_numer = ka**2 * (ka - ki) * (c + ki)
    wt_denom = (c + ka) * ((c + ka)**2 * ki**2 +\
                 np.exp(-ep_ai) * ka**2 * (c + ki)**2)
    mut_numer = ka_mut**2 * (ka_mut - ki_mut) * (c + ki_mut)
    mut_denom = (c + ka_mut) * ((c + ka_mut)**2 * ki_mut**2 +\
                 np.exp(-ep_ai) * ka_mut**2 * (c + ki_mut)**2)
    return 2 * np.exp(-ep_ai) * ((mut_numer / mut_denom) - (wt_numer / wt_denom))

# Define the parameter ranges
theta_range = np.logspace(-5, 5, 500)
theta1 = np.logspace(-3, 3, 500)
theta2 = np.logspace(-3.5, 2.5, 500)
c_range = np.logspace(-5, 5, 500)
c, t = np.meshgrid(c_range, theta_range)
_, t1= np.meshgrid(c_range, theta1)
_, t2 = np.meshgrid(c_range, theta2)
theta_delF = delF(c, t, t)
theta1_delF = delF(c, t1, 1)
theta2_delF = delF(c, 1, t2)

 
theta_deriv = derivF(c, t, t)
theta1_deriv = derivF(c, t1, 1)
theta2_deriv = derivF(c, 1, t2)

#
#%% Instantiate the figure and format the axes
fig, ax = plt.subplots(2, 3, figsize=(6, 4))
for a in ax.ravel():
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.grid(False)

delF1 = ax[0, 0].imshow(theta_delF, origin='bottom', cmap='viridis')
delF2 = ax[0, 1].imshow(theta1_delF, origin='bottom', cmap='viridis')
delF3 = ax[0, 2].imshow(theta2_delF, origin='bottom', cmap='viridis')
ax[0, 0].contour(theta_delF, levels=[-5, -2, 0, 2, 5], colors='w', linestyles='-')
ax[0, 1].contour(theta1_delF, levels=[-5, -2, 0, 2, 5], colors='w', linestyles='-')
ax[0, 2].contour(theta2_delF, levels=[-5, -2, 0, 2, 5], colors='w', linestyles='-')



