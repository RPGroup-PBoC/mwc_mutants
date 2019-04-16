# -*- coding: utf-8 -*- 
import sys
sys.path.insert(0, '../../')
import pandas as pd
import numpy as np
import mut.thermo
import mut.stats
constants = mut.thermo.load_constants()
mut.viz.plotting_style()

# Define the constants relative for drawing samples
n_draws = 800
c_range = np.logspace(-2, 4, 12)
c_range[0] = 0
model_names = ['KaKi_only', 'KaKi_epAI']

# Define the analytical prior distributions
epk_range = np.linspace(-10, 8.52);
epai_range = np.linspace(-10, 10)

ep_a = np.random.normal(0, 5, n_draws)
ep_i = np.random.normal(0, 5, n_draws)
ep_ai = np.random.normal(0, 5, n_draws)
ka = np.exp(ep_a)
ki = np.exp(ep_i)

dfs = []
for m in model_names:

    if m == 'KaKi_only':
        args = dict(ka=ka, ki=ki, ep_ai=constants['ep_AI'])
    else:
        args = dict(ka=ka, ki=ki, ep_ai=ep_ai)

    for i, c in enumerate(c_range):
        arch = mut.thermo.SimpleRepression(R=260, ep_r=-13.9, effector_conc=c, 
                                           **args).fold_change()
        _df = pd.DataFrame([]) #, columns=['IPTGuM', 'ep_a', 'ep_i', 'ka', 'ki', 'ep_ai', 'fold_change'])
        _df['ep_a'] = ep_a
        _df['ka'] = ka
        _df['ep_i'] = ep_i
        _df['ki'] = ki
        _df['ep_ai'] = args['ep_ai']
        _df['fold_change'] = arch
        _df['draw'] = np.arange(n_draws)
        _df['model'] = m
        _df['IPTGuM'] = c_range[i]
        dfs.append(_df)
    
df = pd.concat(dfs)
df.to_csv('../../data/csv/IND_prior_predictive_checks.csv', index=False)




# ##############################################################################
# FIGURE INSTANTIATION 
# ##############################################################################
fig = plt.figure(figsize=(6, 4))
gs = gridspec.GridSpec(4, 6)
ax1 = fig.add_subplot(gs[0:2, 0:2])
ax2 = fig.add_subplot(gs[0:2, 2:4])
ax3 = fig.add_subplot(gs[0:2, 4:])
ax4 = fig.add_subplot(gs[2:, 0:3])
ax5 = fig.add_subplot(gs[2:, 3:])
ax = [ax1, ax2, ax3, ax4, ax5]
for a