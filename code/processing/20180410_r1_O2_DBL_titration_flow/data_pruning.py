# -*- coding: utf-8 -*-
import pandas as pd

# Get the csv file. 
file = pd.read_csv('output/20180410_r1_DBL_fold_change.csv')

# Drop the spurious concentration
file = file[file['mutant'] != 'Y20I-Q294K']

# Resave. 
file.to_csv('output/20180410_r1_DBL_fold_change.csv')