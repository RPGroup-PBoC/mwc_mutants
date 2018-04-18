#! /usr/bin/env python
# -*- coding: utf-8 -*-
import glob
import pandas as pd

# Define the trimmed FC bounds
fc_bounds = (-0.2, 1.3)

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

# Get all valid files in processing folder.
files = glob.glob('processing/2018*/output/*.csv')
all_data = pd.concat([pd.read_csv(f, comment="#") for f in files])

# Drop unnecessary columns.
dropped_cols = [c for c in all_data.keys() if c not in col_names]
all_data.drop(columns=dropped_cols, inplace=True)

# Place the class.
for k in class_dict.keys():
    all_data.loc[all_data['mutant'] == k, 'class'] = class_dict[k]

# Include the method information.
microscopy_files = glob.glob('processing/*microscopy*')
microscopy_dates = [f.split('/')[-1].split('_')[0] for f in microscopy_files]
all_data.loc[:, 'method'] = 'flow cytometry'
for d in microscopy_dates:
    all_data.loc[(all_data['date'] == int(d)), 'method'] = 'microscopy'
    all_data.loc[(all_data['date'] == int(d)), 'IPTGuM'] = 0
    all_data.loc[(all_data['date'] == int(d)), 'operator'] = 'O2'

# Filter the fold change.
filtered_data = all_data[(all_data['fold_change'] >= fc_bounds[0]) &\
                         (all_data['fold_change'] <= fc_bounds[1]) &\
                         (all_data['repressors'] > 0) & (all_data['mutant'] != 'auto') &
                         (all_data['mutant'] != 'delta')]

# Save it to the data directory.
filtered_data.to_csv('../data/csv/compiled_data.csv', index=False)