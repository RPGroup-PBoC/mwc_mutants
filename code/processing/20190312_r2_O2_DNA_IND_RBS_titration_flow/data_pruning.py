#-*- coding: utf-8 -*-
import pandas as pd
file = pd.read_csv('output/20190312_r2_fold_change.csv')
file = file[file['IPTGuM'] != 75]
file.to_csv('output/20190312_r2_fold_change.csv', index=False)
