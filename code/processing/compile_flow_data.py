import glob
import pandas as pd

#============================================================================== 
# list directories with IPTG titration
files = glob.glob('*flow/output/*IPTG*.csv')

# Import them as a single data frame
df = pd.concat(pd.read_csv(f, comment='#') for f in files)

# Export data frame as flow_master.csv
df.to_csv('../../data/flow_master.csv')

# list directories with lacI titration
files = glob.glob('*flow/output/*lacI*.csv')

# Import them as a single data frame
df = pd.concat(pd.read_csv(f, comment='#') for f in files)

# Export data frame as flow_master.csv
df.to_csv('../../data/flow_lacI_master.csv')
