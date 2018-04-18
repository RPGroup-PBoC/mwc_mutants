import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pboc.plotting
import os
import pboc.transcription as trans
import glob
colors = pboc.plotting.set_plotting_style()

# Define the experimental parameters
DATE = 20171205
MUT_CLASS = 'DNA'
OPERATOR = 'O1'


# Load the data.
files = glob.glob(
    'output/{0}_{1}_*microscopy_measurements.csv'.format(DATE, OPERATOR))[0]
data = pd.read_csv(files, comment='#')

# Compute the mean auto and delta.
mean_auto = data[data['rbs'] == 'auto']['mean_intensity'].mean()
mean_delta = data[data['rbs'] == 'delta']['mean_intensity'].mean()

# Group by the mutant
rep_data = data.loc[(data['repressors'] > 0), ['date', 'username', 'IPTG_uM', 'mutant', 'mut_class',
                                               'operator', 'rbs', 'repressors',
                                               'mean_intensity']]
grouped = rep_data.groupby(['mutant', 'repressors'])
dfs = []
for g, d in grouped:
    mean_int = d['mean_intensity'].mean()
    fc = (mean_int - mean_auto) / (mean_delta - mean_auto)
    samp_dict = dict(date=d['date'].unique()[0],
                     username=d['username'].unique()[0],
                     IPTG_uM=d['IPTG_uM'].unique(), mutant=g[0],
                     mut_class=d['mut_class'].unique()[0],
                     operator=d['operator'].unique()[0],
                     rbs=d['rbs'].unique()[0],
                     repressors=g[1], mean_intensity=mean_int,
                     fold_change=fc)
    dfs.append(pd.DataFrame(samp_dict))
fc_df = pd.concat(dfs, axis=0)


# Set the theoretical line for wild-type.
op_en = {'O1': -15.3, 'O2': -13.9, 'O3': -9.7}
ka = 139e-6
ki = 0.53e-6
ep_ai = 4.5
ep_r = op_en[OPERATOR]
effector_conc = 0
R = np.logspace(0, 4, 500)

# Instantiate the architecture.
arch = trans.SimpleRepression(R, ep_r, ka=ka, ki=ki, ep_ai=ep_ai,
                              effector_conc=effector_conc)
theo_fc = arch.fold_change()

# Generate the plot.
fig, ax = plt.subplots(1, 1)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('number of repressors')
ax.set_ylabel('fold-change')

# Plot the theoretical expectation for wt.
_ = ax.plot(R, theo_fc, color='k', label='theory for WT')

# Plot the data
grouped = fc_df.groupby('mutant')
for g, d in grouped:
    _ = ax.plot(d['repressors'], d['fold_change'], 'o', label=g)

_ = ax.legend(loc='lower left')


plt.tight_layout()
plt.savefig('output/{0}_{1}_{2}_repressor_titration.png'.format(DATE, MUT_CLASS,
                                                                OPERATOR),
            bbox_inches='tight')

target = 'output/{0}_{1}_{2}_repressor_titration_microscopy_foldchange.csv'.format(
    DATE, MUT_CLASS, OPERATOR)
fc_df.to_csv('output/tmp.csv')
with open('output/tmp.csv', 'r+') as f:
    csv = f.readlines()
with open('comments.txt', 'r+') as f:
    header = f.readlines()
with open(target, 'w') as f:
    for line in header:
        f.write(line)
    for line in csv:
        f.write(line)
os.remove('output/tmp.csv')
