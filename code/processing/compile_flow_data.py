import glob
import pandas as pd

#============================================================================== 
# list directories
files = glob.glob('*flow/output/*.csv')

# Import them as a single data frame
df = pd.concat(pd.read_csv(f, comment='#') for f in files)

# Export data frame as flow_master.csv
df.to_csv('../../data/flow_master.csv')
