# -*- coding: utf-8 -*-
"""
This is a representative script for processing raw flow cytometry files. It is
generic to all processing steps and requires only minor changes to parameters 
at the front of the script to run on any data set. 
"""
import numpy as np
import pandas as pd
import glob
import imp
import mut.flow

# Define the experiment parameters
DATE = 20190410 # Date of the experiment
RUN_NO = 2 # Integer identifier of the run number for the day
USERNAME = 'gchure' # Identifier of who took the data
gating_fraction = 0.4 # Sets the percentile of the flow data to keep.

# Load all files. 
files = glob.glob('../data/flow/{0}*_r{1}*.csv'.format(DATE, RUN_NO))

# Set up the DataFrame
colnames = ['date', 'username', 'mutant', 'operator', 'strain', 'IPTGuM',
            'mean_FITC_H']
df = pd.DataFrame([], columns=colnames)

for f in files:
    # Get the identifying finformation.
    date, _, operator, strain, mutant, conc = f.split('/')[-1].split('_')
    conc = float(conc.split('uM')[0])
    rep = int(strain.split('R')[-1])

    # Load in the data
    data = pd.read_csv(f)
    gated = mut.flow.gaussian_gate(data, gating_fraction)

    # Compute the mean
    mean_FITC = gated['FITC-H'].mean()

    # Assemble the dictionary
    samp_dict = dict(date=date, username=USERNAME, mutant=mutant,
                     operator=operator, strain=strain, IPTGuM=conc,
                     mean_FITC_H=mean_FITC, repressors=rep)
    df = df.append(samp_dict, ignore_index=True)


fc_dfs = []
grouped = df[df['mutant']!='auto'].groupby(['IPTGuM', 'operator'])
mean_auto_df = df[df['mutant'] == 'auto']
for g, d in grouped:
    mean_auto = mean_auto_df[mean_auto_df['IPTGuM']==g[0]]['mean_FITC_H'].values[0]
    mean_delta = d.loc[d['mutant'] == 'delta']['mean_FITC_H'].values[0]
    d['fold_change'] = (d['mean_FITC_H'] - mean_auto) /  (mean_delta - mean_auto)
    fc_dfs.append(d)

fold_change_df = pd.concat(fc_dfs, axis=0)

# Save to a CSV.
fold_change_df.to_csv(
    './output/{0}_r{1}_fold_change.csv'.format(DATE, RUN_NO))