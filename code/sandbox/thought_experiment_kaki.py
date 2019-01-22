# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
colors = mut.viz.personal_style()
constants = mut.thermo.load_constants()

# Define the reference state. 
c_0 = 50
R_0 = 260
ref_state = mut.thermo.SimpleRepression(R=R_0, ep_r=constants['O2'], 
                                  effector_conc=c_0, ka=constants['Ka'],
                                  ki=constants['Ki'], n_sites=constants['n_sites'],
                                  ep_ai=constants['ep_AI'])

c_range = np.logspace(-2, 4, 500)
ref_strain = mut.thermo.SimpleRepression(R=R_0, ep_r=constants['O2'],
                                        effector_conc=c_range, ka=constants['Ka'],
                                        ki=constants['Ki'], ep_ai=constants['ep_AI'],
                                        n_sites=constants['n_sites'])


ka_range = np.logspace(-3, 3, 500)
ki_range = np.logspace(-3, 3, 500)
ka, ki, c = np.meshgrid(ka_range, ki_range, c_range)

# Generate the fold-change data set. 
exp_arch = mut.thermo.SimpleRepression(R=R_0, ep_r=constants['O2'],
                                      effector_conc=c, ka=ka, 
                                      ki=ki, n_sites=constants['n_sites'],
                                      ep_ai=constants['ep_AI'])

exp_fc = exp_arch.fold_change()

# Compute the empirical delta F. 
measured_F = np.log((1 / exp_fc) - 1)

# Compute the empDeltaF
emp_dF = -ref_state.bohr_parameter() - measured_F
 # Compute the ddf 
ddf = emp_dF - (-ref_strain.bohr_parameter() + ref_state.bohr_parameter())

# Look at the evolution using holoviews. 
i = 100 
plt.imshow(ddf[i, :,  :],  cmap='bone')
c = plt.contour(ddf[i, :,  :], levels=[-20, -10, -5, 0, 5, 10, 15, 20], colors='w', linestyles='-')
plt.clabel(c)
# plt.colorbar()


# plt.plot(np.arange(500), ddf[-20, 10, :])
