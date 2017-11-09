import numpy as np
import matplotlib.pyplot as plt
import pboc.plotting
import pandas as pd
import glob
colors = pboc.plotting.set_plotting_style(return_colors=True)

# Set the experiment constants.
DATE = 20171107
RUN_NO = 1
MUTANT = 'antilac'

# Load the data set.
fc_file = glob.glob('output/*fold_change.csv')[0]
data = pd.read_csv(fc_file)

# Instantiate the figure canvas
fig, ax = plt.subplots(1, 1, figsize=(6, 4))

# Add labels and scaling
ax.set_xlabel('IPTG [M]')
ax.set_ylabel('fold-change')
ax.set_xscale('log')

# Group the data by operator
grouped = data.loc[data['strain'] == 'RBS1027'].groupby('operator')

# Plot the inensity curves.
color_id = 0
for g, d in grouped:
    _ = ax.plot(d['IPTGuM'] / 1E6, d['fold_change'], '--o',
                color=colors[color_id], label=g)
    color_id += 1

# Add a legend.
_ = ax.legend(loc='upper right', title='operator')

# Save the figure.
plt.savefig('output/{0}_r{1}_{2}_fold_change_curve.png'.format(DATE, RUN_NO,
                                                               MUTANT))
