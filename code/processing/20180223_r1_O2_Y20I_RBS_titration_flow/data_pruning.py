# -*- coding: utf-8 -*-
import pandas as pd

# Get the csv file. 
file = pd.read_csv('output/20180223_r1_DNA_fold_change.csv')

# Drop the spurious concentration
file = file[file['IPTGuM'] != 50]

# Resave. 
file.to_csv('output/20180223_r1_DNA_fold_change.csv')