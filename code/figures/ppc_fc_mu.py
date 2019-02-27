#-*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import mut.stats
import mut.viz
mut.viz.plotting_style()
pboc = mut.viz.color_selector('pboc')


# Load the observed data
data = pd.read_csv('../../data/csv/compiled_data.csv')
data.dropna(inplace=True)
data.groupby(['mutant', 'IPTGuM', 'operator', 'repressors']).count()


# Load the posterior predictive checks. 
ppc = pd.read_csv('../../data/csv/inferred_fc_ppc.csv')

fig, ax = plt.subplots(1, 1)

# Plot the credible region. 
y = np.linspace(0, 1, len(ppc))
x_min = np.sort(ppc['hpd_min'])
x_max = np.sort(ppc['hpd_max'])
ax.fill_betweenx(y, x_min, x_max, color=pboc['light_red'], step='pre')

# Plot the data. 
x_data = np.sort(data['fold_change'])
ax.step(x_data, y, 'k', lw=0.75)



v = data[(data['operator']=='O3') & (data['mutant']=='Q294V')]

fig, ax = plt.subplots()
ax.semilogx(v['IPTGuM'], v['fold_change'], 'o')
ax.set_ylim([-0.15, 1.15])