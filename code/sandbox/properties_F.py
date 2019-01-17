# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.viz 
import mut.thermo
colors = mut.viz.personal_style()
constants = mut.thermo.load_constants()

rep_range = np.logspace(0, 4, 200)
dr = (1 + np.exp(-constants['ep_AI'])) * (rep_range/4.6E6) * np.exp(13.9)

mwc_arch = mut.thermo.MWC(effector_conc=10, ka=constants['Ka'],
                      ki=constants['Ki'], ep_ai=constants['ep_AI'],
                      n_sites=constants['n_sites'])

F_sat = np.log(mwc_arch.saturation())
F_leak = np.log(mwc_arch.leakiness())


F_sat = np.linspace(-10, 10, 100)
F_leak = np.linspace(-10, 10, 100)

sat, leak = np.meshgrid(F_sat, F_leak)

c1 = plt.matshow(sat - leak)
cl1 = plt.contour(sat - leak, levels=[-10, -5, 0, 5, 10], colors='w',
                  linestyles='-')
plt.clabel(cl1)

cl2 = plt.contour((1 + np.exp(sat))**-1 - (1 + np.exp(leak))**-1,
                 levels=[-0.1, 0.0, 0.01, 0.1, 0.5, 1], colors='w', linestyles='-')


c2 = plt.matshow(np.log(np.abs((1 + np.exp(sat))**-1 - (1 + np.exp(leak))**-1)),
                cmap='viridis')
cl1 = plt.contour(sat - leak, levels=[-10, -5, 0, 5, 10], colors='w',
                  linestyles='-')
plt.clabel(cl1)