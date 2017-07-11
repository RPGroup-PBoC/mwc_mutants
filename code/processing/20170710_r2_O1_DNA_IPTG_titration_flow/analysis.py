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
date = 20170710
run = 'r2'
operator = 'O1'
mutant_class = 'dna' # options dna/inducer/double

# Read the CSV file with the mean fold change
df = pd.read_csv('output/' + str(date) + '_' + run + '_' + operator + '_' + \
          mutant_class + '_IPTG_titration_MACSQuant.csv', comment='#')
rbs = df.strain.unique()

#==============================================================================
# Plot all raw data
#==============================================================================

plt.figure()
for strain in rbs[np.array([r != 'auto' and r != 'delta' for r in rbs])]:
    plt.plot(df[df.strain == strain].sort_values(by='IPTG_uM').IPTG_uM * 1E-6,
             df[df.strain == strain].sort_values(by='IPTG_uM').fold_change,
             marker='o', linewidth=1, linestyle='--', label=strain)
plt.xscale('log')
plt.xlabel('IPTG (M)')
plt.ylabel('fold-change')
plt.ylim([-0.01, 1.2])
plt.xlim([1E-8, 1E-2])
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig('output/' + str(date) + '_' + operator + '_' + mutant_class + \
            '_IPTG_titration_data.png')

#==============================================================================
# Plot the WT control with the parameters determined from the MWC induction
# project
#==============================================================================
# Define the IPTG concentrations to evaluate
IPTG = np.logspace(-7, -2, 100)
IPTG_lin = np.array([0, 1E-7])

energies = {'O1': -15.3, 'O2': -13.9, 'O3': -9.7, 'Oid': -17}
ka = -np.log(139.59)
ki = -np.log(0.53)

df_wt = df[df.strain=='wt']

# Initialize figure
plt.figure()
# Plot theoretical prediction
# Log scale
plt.plot(IPTG, mwc.fold_change_log(IPTG * 1E6,
    ka=ka, ki=ki, epsilon=4.5,
    R=df_wt.repressors.unique()[0],
    epsilon_RA=energies[df_wt.operator.unique()[0]]),
    color='black', label='prediction')
# Linear scale
plt.plot(IPTG_lin, mwc.fold_change_log(IPTG_lin * 1E6,
    ka=ka, ki=ki, epsilon=4.5,
    R=df_wt.repressors.unique()[0],
    epsilon_RA=energies[df_wt.operator.unique()[0]]),
    linestyle='--', color='black', label=None)

# Plot data
plt.scatter(df_wt.IPTG_uM / 1E6, df_wt.fold_change, label='wt data')

plt.xscale('symlog', linthreshx=1E-7, linscalex=0.5)
plt.xlabel('IPTG (M)', fontsize=15)
plt.ylabel('fold-change', fontsize=16)
plt.ylim([-0.05, 1.1])
plt.xlim([-5E-9, np.max(IPTG)])
plt.legend(loc='upper left')
plt.tick_params(labelsize=14)
plt.tight_layout()
plt.savefig('output/' + str(date) + '_' + operator + '_' + mutant_class + \
            '_wt_titration.png')
