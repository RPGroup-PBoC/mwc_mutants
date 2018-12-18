import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import sys
sys.path.insert(0, '../../../')
import mut.viz
import mut.thermo
mut.viz.plotting_style()
colors =list(mut.viz.color_selector('pboc').values())[::1]

# Set the experiment constants.
DATE = 20181206
RUN_NO = 1
MUTANT = 'IND'

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
grouped = fc.groupby(['mutant', 'repressors', 'operator'])

# Plot the inensity curves.
color_id = 0
for g, d in grouped:
    _ = ax.plot(d['IPTGuM'] / 1E6, d['fold_change'], '--o',
                color=colors[color_id], label=g)
    color_id += 1

# Add a legend.
_ = ax.legend(loc='upper left', title='operator')
ax.set_ylim([-0.1, 1.5])
# Save the figure.
plt.savefig('output/{0}_r{1}_{2}_fold_change_curve.png'.format(DATE, RUN_NO,
                                                               MUTANT))
