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
             'operator', 'repressors', 'mutant', 'fold_change']

# Generate list of "Accepted" experiments.
accepted = []
files = glob.glob('processing/2018*')
for _, d in enumerate(files):
    info = mut.io.scrape_frontmatter(f'{d}')
    if info['status'].lower() == 'accepted':
        accepted.append(glob.glob(f'{d}/output/*.csv')[0])

# Get all valid files in processing folder.
all_data = pd.concat([pd.read_csv(f, comment="#") for f in accepted],sort=False)

# Drop unnecessary columns.
dropped_cols = [c for c in all_data.keys() if c not in col_names]
all_data.drop(columns=dropped_cols, inplace=True)

# Place the class.
for k in class_dict.keys():
    all_data.loc[all_data['mutant'] == k, 'class'] = class_dict[k]

# Include the method information.
# microscopy_files = glob.glob('processing/*microscopy*')
# microscopy_dates = [f.split('/')[-1].split('_')[0] for f in microscopy_files]
# all_data.loc[:, 'method'] = 'flow cytometry'
# for d in microscopy_dates:
#     all_data.loc[(all_data['date'] == int(d)), 'method'] = 'microscopy'
# #     all_data.loc[(all_data['date'] == int(d)), 'IPTGuM'] = 0
#     all_data.loc[(all_data['date'] == int(d)), 'operator'] = 'O2'

# Filter the fold change.
filtered_data = all_data[(all_data['repressors'] > 0) & (all_data['mutant'] != 'auto') &
                         (all_data['mutant'] != 'delta')]
                         
# Save it to the data directory.
filtered_data.to_csv('../data/csv/compiled_data.csv', index=False)

# Compute the summarized statistics and save. 
summarized = filtered_data.groupby(['IPTGuM', 'class', 'operator', 'repressors', 'mutant']).apply(mut.stats.compute_mean_sem)
summarized = pd.DataFrame(summarized).reset_index()
summarized.to_csv('../data/csv/summarized_data.csv', index=False)