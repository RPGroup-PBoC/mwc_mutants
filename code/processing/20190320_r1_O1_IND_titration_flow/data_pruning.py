import pandas as pd
file = pd.read_csv('output/20190320_r1_fold_change.csv')
file = file[file['IPTGuM'] != 75]
file.to_csv('output/20190320_r1_fold_change.csv', index=False)