# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pystan
sys.path.insert(0, '../../')
import mut.viz
import mut.bayes
import mut.stats
import mut.thermo
colors = mut.viz.pub_style()

# Load the data.
data = pd.read_csv('../../data/csv/compiled_data.csv')
wt_data = data[data['mutant'] == 'wt']

# Load the stan model for the global fit.
model_code = mut.bayes.assemble_StanModelCode(
    '../stan/hierarchical_global_fit.stan', '../stan/functions.stan')
model = pystan.StanModel(model_code=model_code)

# %% Assemble the data dictionary and sample.
data_dict = {'J': 1, 'N': len(wt_data), 'trial': np.ones(len(wt_data)).astype(int),
             'R': wt_data['repressors'], 'n_ns': 4.6E6, 'ep_AI': 4.5, 'n_sites': 2,
             'c': wt_data['IPTGuM'], 'fc': wt_data['fold_change']}
samples = model.sampling(data=data_dict, iter=10000, chains=4)

# %% Process the sampling data.
samples_df = mut.bayes.chains_to_dataframe(samples)
sample_stats = mut.stats.compute_statistics(samples_df)

# Save to csv
samples_df.to_csv('../../data/csv/WT_global_fit_samples.csv', index=False)
sample_stats.to_csv('../../data/csv/WT_global_fit_parameters.csv', index=False)

# %% Plot the compiled WT data.
fig, ax = plt.subplots(1, 1, figsize=(6, 4))
ax.set_xlabel('IPTG [µM]', fontsize=10)
ax.set_ylabel('fold-change', fontsize=10)
ax.set_xscale('log')
ax.set_xlim([1E-8, 1E-2])

# Plot the data.
grouped = wt_data.groupby('IPTGuM').apply(mut.stats.compute_mean_sem)
fc_mean = pd.DataFrame(grouped).reset_index()
_ = ax.errorbar(fc_mean['IPTGuM']/1E6, fc_mean['mean'], fc_mean['sem'], fmt='o', color=colors['red'],
                ms=4, lw=1, markerfacecolor='w', markeredgecolor=colors['red'], markeredgewidth=1)

# Extract the parameter values
modes = {}
for i, p in enumerate(sample_stats['parameter'].unique()):
    if '_tilde' not in p:
        modes[p] = sample_stats[sample_stats['parameter'] == p]['mode'].values[0]
c_range = np.logspace(-8, -2, 500)

# Compute the best-fit line
arch = mut.thermo.SimpleRepression(
    R=260, ep_r=modes['ep_R'], effector_conc=c_range, ep_ai=4.5, ka=modes['ka'] / 1E6, ki=modes['ki'] / 1E6)
fc_theo = arch.fold_change()

# Compute the credible region.
cred_region = np.zeros((2, len(c_range)))
for i, c in enumerate(c_range):
    arch = mut.thermo.SimpleRepression(
        R=260, ep_r=samples_df['ep_R'], effector_conc=c, ep_ai=4.5, ka=samples_df['ka'] / 1E6, ki=samples_df['ki'] / 1E6)
    prob = arch.fold_change()
    cred_region[:, i] = mut.stats.compute_hpd(prob, 0.95)

_ = ax.plot(c_range, fc_theo, color=colors['red'], lw=1)
_ = ax.fill_between(
    c_range, cred_region[0, :], cred_region[1, :], color=colors['light_red'], alpha=0.75)
plt.tight_layout()
plt.savefig('../../figures/wt_global_fit.pdf', bbox_inches='tight', dpi=300)
