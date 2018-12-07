# -*- coding: utf-8 -*-
import pandas as pd

# Get the csv file. 
file = pd.read_csv('output/20180522_r1_IND_fold_change.csv')

# Drop the spurious concentration
file = file[(file['run_no'] != 'r5') & (file['IPTGuM'] != 100)]

# Resave. 
file.to_csv('output/20180522_r1_IND_fold_change.csv')