import pandas as pd
file = pd.read_csv('output/20190319_r2_fold_change.csv')
file = file[file['IPTGuM'] < 1000]
file.to_csv('output/20190319_r2_fold_change.csv', index=False)