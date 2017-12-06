import glob
import pandas as pd

#==============================================================================
# # list directories with IPTG titration
# files = glob.glob('*microscopy/output/*IPTG*.csv')
#
# # Import them as a single data frame
# df = pd.concat(pd.read_csv(f, comment='#') for f in files)
# # Rename the column class to mut_class to avoid conflicts
# df.rename(columns={'class': 'mut_class'}, inplace=True)
#
# # Export data frame as flow_master.csv
# df.to_csv('../../data/flow_master.csv')

# list directories with lacI titration
print('compiling repressor titration files...')
files = \
glob.glob('*repressor_titration_microscopy/output/*microscopy_foldchange.csv')

# Import them as a single data frame
df = pd.concat(pd.read_csv(f, comment='#') for f in files)

# Export data frame as
df.to_csv('../../data/microscopy_repressor_master.csv')

print('done!')
