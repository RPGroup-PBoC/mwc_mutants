import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import sys
sys.path.insert(0, '../../../')
import mut.viz
import mut.thermo
colors = mut.viz.pub_style(return_colors=True)
colors = list(colors.values())[::3]

# Set the experiment constants.
DATE = 20180321
RUN_NO = 3
MUTANT = 'DNA'
OPERATOR = 'O2'
op_en = {'O1': -15.3, 'O2': -13.9, 'O3': -9.7}

# Load the data set.
fc_file = glob.glob('output/*fold_change.csv')[0]
data = pd.read_csv(fc_file)

# Instantiate the figure canvas
fig, ax = plt.subplots(1, 1)

# Add labels and scaling
ax.set_xlabel('IPTG [M]')
ax.set_ylabel('fold-change')
ax.set_xscale('log')

# Group the data by operator
# Remove auto and delta.
fc = data.loc[(data['mutant'] != 'auto') & (data['mutant'] != 'delta')]
grouped = fc.groupby(['mutant', 'repressors'])

# Plot the inensity curves.
color_id = 0
for g, d in grouped:
    _ = ax.plot(d['IPTGuM'] / 1E6, d['fold_change'], '--o',
                color=colors[color_id], label=g)
    color_id += 1

# Add a legend.
_ = ax.legend(loc='upper left', title='operator')

ax.set_ylim([-0.1, 1.2])
# Save the figure.
plt.savefig('output/{0}_r{1}_{2}_fold_change_curve.png'.format(DATE, RUN_NO,
                                                               MUTANT))
# Compute the expected curve from WT.
ka = 139e-6
ki = 0.53e-6
ep_AI = 4.5
ep_RA = op_en[OPERATOR]
c_range = np.logspace(-8, -2, 500)
R = 260

# Instantiate the architecture and compute the fold-change.
arch = mut.thermo.SimpleRepression(R, ep_RA, ka=ka, ki=ki,
                                   ep_ai=ep_AI, effector_conc=c_range)
fc_theo = arch.fold_change()

# Generate the figure.
fig, ax = plt.subplots(1, 1)

ax.set_xlabel('IPTG [M]')
ax.set_ylabel('fold-change')
ax.set_xscale('log')

# Plot the data.
wt_data = fc.loc[fc['mutant'] == 'wt']
_ = ax.plot(wt_data['IPTGuM'] / 1E6,
            wt_data['fold_change'], 'o', label='WT data')

# Plot the prediction
_ = ax.plot(c_range, fc_theo, 'k-', label='WT prediction')
plt.legend(loc='upper left')
ax.ylim([-0.1, 1.2])
plt.tight_layout()
plt.savefig(
    'output/{0}_r{1}_{2}_WT_titration.png'.format(DATE, RUN_NO, MUTANT))
