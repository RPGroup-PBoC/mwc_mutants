# Import dependencies
import numpy as np
import pandas as pd

# Import matplotlib stuff for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Seaborn, useful for graphics
import seaborn as sns

# Set the plotting style.
import sys
sys.path.insert(0, '../../')

import mwc_mutants_utils as mwc
mwc.set_plotting_style()


# Define variables to use over the script
date = 20170616
run = 'r1'
operator = 'O1'
mutant_class = 'dna' # options dna/inducer/double

# Read the CSV file with the mean fold change
df = pd.read_csv('output/' + str(date) + '_' + run + '_' + operator + '_' + \
          mutant_class + '_lacI_titration_MACSQuant.csv', comment='#')
rbs = df.strain.unique()

#==============================================================================
# Plot all raw data
#==============================================================================

plt.figure()
for strain in rbs[np.array([r != 'auto' and r != 'delta' for r in rbs])]:
    plt.plot(df[df.strain == strain].sort_values(by='repressors').repressors,
             df[df.strain == strain].sort_values(by='repressors').fold_change,
             marker='o', linewidth=1, linestyle='--', label=strain)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('repressors/cell')
plt.ylabel('fold-change')
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig('output/' + str(date) + '_' + operator + '_' + mutant_class + \
            '_lacI_titration_data.png')
