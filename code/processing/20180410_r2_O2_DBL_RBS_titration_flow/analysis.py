import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import sys
sys.path.insert(0, '../../../')
import mut.viz
import mut.thermo
colors = mut.viz.pub_style(return_colors=True)

colors = list(colors.values())[::2]

# Set the experiment constants.
DATE = 20180410
RUN_NO = 2
MUTANT = 'DBL'
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

# Save the figure.
plt.savefig('output/{0}_r{1}_{2}_fold_change_curve.png'.format(DATE, RUN_NO,
                                                               MUTANT))
