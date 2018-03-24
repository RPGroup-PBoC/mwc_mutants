import numpy as np
import pandas as pd
import pboc.flow
import glob
import imp
imp.reload(pboc.flow)
# Define the experiment parameters
DATE = 20171107
RUN_NO = 1
OPERATORS = ['O1', 'O2', 'O3']
USERNAME = 'gchure'
mutant = 'antilac'
gating_fraction = 0.4

# Load all files.
files = glob.glob('../../../data/flow/csv/{0}_r{1}*.csv'.format(DATE, RUN_NO))
np.sort(files)
# Set up the DataFrame
colnames = ['date', 'username', 'mutant', 'operator', 'strain', 'IPTGuM',
            'mean_FITC_H']
df = pd.DataFrame([], columns=colnames)
for f in files:
    f
    # Get the identifying finformation.
    try:
        date, _, mutant, operator, strain, conc = f.split('/')[-1].split('_')
    except:
        date, _, mutant, strain, conc = f.split('/')[-1].split('_')
        operator = None
    conc = float(conc.split('uM')[0])

    # Load in the data
    data = pd.read_csv(f)
    gated = pboc.flow.gaussian_gate(data, gating_fraction)

    # Compute the mean
    mean_FITC = gated['FITC-H'].mean()

    # Assemble the dictionary
    samp_dict = dict(date=date, username=USERNAME, mutant=mutant,
                     operator=operator, strain=strain, IPTGuM=conc,
                     mean_FITC_H=mean_FITC)
    df = df.append(samp_dict, ignore_index=True)

# Compute the fold-change
mean_auto = df.loc[df['strain'] == 'auto']
auto_grouped = mean_auto.groupby('IPTGuM')
auto_dict = {}
for g, d in auto_grouped:
    auto_dict[g] = d['mean_FITC_H'].values[0]

fc_dfs = []
grouped = df.groupby(['operator', 'IPTGuM'])
for g, d in grouped:
    d
    mean_delta = d.loc[d['strain'] == 'delta']['mean_FITC_H'].values
    d['fold_change'] = (d['mean_FITC_H'] - auto_dict[g[1]]) / \
        (mean_delta - auto_dict[g[1]])
    fc_dfs.append(d)

fold_change_df = pd.concat(fc_dfs, axis=0)

# Save to a CSV.
fold_change_df.to_csv(
    'output/{0}_r{1}_{2}_fold_change.csv'.format(DATE, RUN_NO, mutant))


# Add the comments and save to the data/csv file.
target = '../../../data/csv/{0}_r{1}_{2}_fold_change.csv'.format(DATE, RUN_NO,
                                                                 mutant)
with open('comments.txt', 'r') as f:
    comments = f.read().splitlines()

with open(target, 'a') as f:
    for line in comments:
        f.write(line)
    fold_change_df.to_csv(f, index=False)
