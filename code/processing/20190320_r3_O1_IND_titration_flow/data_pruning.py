import pandas as pd
file = pd.read_csv('output/20190320_r3_fold_change.csv')
file.loc[file['repressors'] > 0, 'repressors'] = 260
file.to_csv('output/20190320_r3_fold_change.csv', index=False)