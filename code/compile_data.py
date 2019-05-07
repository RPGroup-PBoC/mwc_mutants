#! /usr/bin/env python
# -*- coding: utf-8 -*-
import glob
import sys
import pandas as pd
sys.path.insert(0, '../')
import mut.io
import mut.stats

# Assemble the idenfier dictionary.
mutants = ['Y20I', 'Q21A', 'Q21M', 'Q294K', 'Q294V', 'Q294R', 'F164T']
class_dict = {'Y20I': 'DNA', 'Q21A': 'DNA', 'Q21M': 'DNA',
              'Q294K': 'IND', 'Q294V': 'IND', 'Q294R': 'IND', 'F164T': 'IND',
              'wt': 'WT'}
keys = class_dict.keys()
for k in mutants:
    for j in mutants:
        if (k != j) and (class_dict[k] == 'DNA') & (class_dict[j] == 'IND'):
            class_dict['{}-{}'.format(k, j)] = 'DBL'

# Define the necessary columns.
col_names = ['date', 'username', 'class', 'IPTGuM',
             'operator', 'repressors', 'mutant', 'fold_change', 'mean_FITC']

# Generate list of "Accepted" experiments.
accepted = []
info_df = pd.DataFrame([])
files = glob.glob('processing/*flow*')
for _, d in enumerate(files):
    info = mut.io.scrape_frontmatter(f'{d}')
    out = d.split('/')[1].split('_')
    run = int(out[1][1])
    info['date'] = out[0]
    info['run'] = run
    if info['status'].lower() == 'accepted':
        accepted.append(glob.glob(f'{d}/output/*.csv')[0])
        info_df = info_df.append(info, ignore_index=True)

info_df = info_df[['date', 'run', 'status']]
info_df.to_csv('../data/csv/valid_data_idx.csv')

# Get all valid files in processing folder.
all_data = pd.concat([pd.read_csv(f, comment="#") for f in accepted],sort=False)

# Drop unnecessary columns.
dropped_cols = [c for c in all_data.keys() if c not in col_names]
all_data.drop(columns=dropped_cols, inplace=True)

# Place the class.
for k in class_dict.keys():
    all_data.loc[all_data['mutant'] == k, 'class'] = class_dict[k]

# Filter the fold change.
filtered_data = all_data[(all_data['repressors'] > 0) & (all_data['mutant'] != 'auto') &
                         (all_data['mutant'] != 'delta')]
                         
# Save it to the data directory.
filtered_data.to_csv('../data/csv/compiled_data.csv', index=False)

# Compute the summarized statistics and save. 
summarized = filtered_data.groupby(['IPTGuM', 'class', 'operator', 'repressors', 'mutant']).apply(mut.stats.compute_mean_sem)
summarized = pd.DataFrame(summarized).reset_index()
summarized.to_csv('../data/csv/summarized_data.csv', index=False)