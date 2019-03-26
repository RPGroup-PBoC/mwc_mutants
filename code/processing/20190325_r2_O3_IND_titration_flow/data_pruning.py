import pandas as pd
file = pd.read_csv('output/20190325_r2_fold_change.csv')
file = file[(file['mutant'] != 'Q294V') & (file['IPTGuM'] != 5000)]
file.to_csv('output/20190325_r2_fold_change.csv')