# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import mut.viz
import mut.thermo
import tqdm
constants = mut.thermo.load_constants()
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load the samples for the DNA binding energy
samples = pd.read_csv('../../data/csv/DNA_binding_energy_samples.csv')
samples = samples[(samples['mutant']=='Q21M') & (samples['repressors']==260)]

# Sort by the logprob and assign a color sequence
cmap = sns.color_palette('viridis', n_colors=len(samples))
samples.sort_values(by='lp__', inplace=True)
samples['color'] = cmap

# Load the fold-change data
data = pd.read_csv('../../data/csv/compiled_data.csv')
data = data[(data['mutant']=='Q21M') & (data['repressors']==260)]

# ##############################################################################
# POSTERIOR PREDICTIVE CHECKS
# ##############################################################################
dfs = [] 
counts = data.groupby('IPTGuM').count()['fold_change'].values.astype(int)
IPTGuM = data['IPTGuM'].unique()
c, ep_r = np.meshgrid(IPTGuM, samples['ep_RA'])
fc_mu = mut.thermo.SimpleRepression(R=260, ep_r=ep_r,
                                    ka=constants['Ka'], ki=constants['Ki'],
                                    effector_conc=c,
                                    ep_ai=constants['ep_AI']).fold_change()

sigma = samples['sigma'].values
for i in tqdm.tqdm(range(len(samples))):
      for j in range(len(IPTGuM)):
        draws = np.random.normal(fc_mu[i,j], sigma[i], size=counts[j])
        _df = pd.DataFrame([], columns = ['IPTGuM'])
        _df['fc_mu'] = draws
        _df['IPTGuM'] = IPTGuM[j]
        _df['samp'] = i
        dfs.append(_df)
ppc_df = pd.concat(dfs)


# ##############################################################################
# FIGURE INSTANTIATION
# ##############################################################################
fig = plt.figure(figsize=(6,  3))
gs = gridspec.GridSpec(4, 9)
ax1 = fig.add_subplot(gs[0:2, 0:2])
ax2 = fig.add_subplot(gs[2:4, 0:2])
ax3 = fig.add_subplot(gs[2:4, 2:4])
ax4 = fig.add_subplot(gs[:, 5:])
ax = [ax1, ax2, ax3, ax4]
for a in ax:
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)

# Turn off axes where not needed
ax1.set_xticklabels([])
ax1.set_yticks([])
ax3.set_yticks([])

# Add labels
ax[1].set_xlabel(r'$\Delta\varepsilon_{RA}$ [$k_BT$]', fontsize=8)
ax[1].set_ylabel('$\sigma$', fontsize=8)
ax[2].set_xlabel('$\sigma$', fontsize=8)
ax[-1].set_xlabel('IPTG [µM]', fontsize=8)
ax[-1].set_ylabel('fold-change', fontsize=8)

# Add appropriate formatting
ax[3].set_xscale('symlog', linxthresh=1E-3)

# Adjust limits
ax[3].set_ylim([-0.15, 0.9])

# Add panel labels
fig.text(0.05, 0.90, '(A)', fontsize=8)
fig.text(0.5, 0.90, '(B)', fontsize=8)

# ##############################################################################
# JOINT DISTRIBUTION
# ##############################################################################
ax[1].scatter(samples['ep_RA'].values, samples['sigma'].values,  c=cmap, marker='.', s=0.5)

# ##############################################################################
# MARGINAL DISTRIBUTIONS
# ##############################################################################
ax[0].hist(samples['ep_RA'].values, bins=30, color=cmap[10], histtype='stepfilled', alpha=0.5)
ax[2].hist(samples['sigma'].values, bins=30, color=cmap[10], histtype='stepfilled', alpha=0.5)

# ##############################################################################
# DATA POINTS
# ##############################################################################
for g, d in data.groupby(['IPTGuM']):
    if g==0:
        label = 'data'
    else:
        label = '__nolegend__'
    ax[3].plot(d['IPTGuM'], d['fold_change'], '.', color='k', ms=5, 
                markerfacecolor='w', alpha=0.5, label=label)

# ##############################################################################
# POSTERIOR PREDICTIVE CHECKS
# ##############################################################################
perc_low = []
perc_high = []
for g, d in ppc_df.groupby('IPTGuM'):
    _low, _high = np.percentile(d['fc_mu'].values, (0.5, 99.5))
    perc_low.append(_low)
    perc_high.append(_high)
ax[3].fill_between(IPTGuM, perc_low, perc_high, color='rebeccapurple', alpha=0.4,
                  label =r'$y \sim \mathcal{N}(\mu, \sigma)$')
ax[3].legend(fontsize=8)
plt.savefig('../../figures/FigSX_epRA_post_pc.pdf', bbox_inches='tight')


