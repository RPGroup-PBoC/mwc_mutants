# -*- coding: utf-8 -*-
import pandas as pd

# Get the csv file. 
file = pd.read_csv('output/20180611_r1_IND_fold_change.csv')

# Drop the spurious concentration
file = file[(file['mutant'] != 'Q294V') & (file['mutant'] !='F164T')]

# Resave. 
file.to_csv('output/20180611_r1_IND_fold_change.csv')