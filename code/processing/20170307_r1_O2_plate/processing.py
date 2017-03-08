# Import dependencies.
import os
import sys
import itertools
import glob
import re
import numpy as np
import pandas as pd
import string

# Set the plotting style.
import sys
sys.path.insert(0, '../../')

# Import matplotlib stuff for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Seaborn, useful for graphics
import seaborn as sns

import mwc_mutants_utils as mwc
mwc.set_plotting_style()

#============================================================================== 
# Variables to edit on each script
#============================================================================== 
# Define variables to use over the script
date = 20170307
username = 'mrazomej'
run = 'r1'
operator = 'O1'
#============================================================================== 
# List the target directory.
datadir = '../../../data/plate_reader/'

# read data
files = np.array(os.listdir(datadir))
# Find files that contain the date, run number and are CSV
csv_bool = np.array([str(date) in f and run in f \
                     and 'xls' in f for f in files])
filename = files[np.array(csv_bool)][0]

df = pd.read_excel(datadir + filename)

#============================================================================== 
# Extract media blank data from the last column
idx_blank = [s for s in df.index if '12' in s]

# Substract mean values from all data
df[['YFP_blank', 'OD_blank']] = df - df.ix[idx_blank].mean()
# Remove blank measurements
df = df.drop(idx_blank)

# make data frame tidy by including repressor, operator and IPTG concentration
df['operator'] = ['O2'] * len(df)

# add strain column
# list strains in order
strains = ['auto', 'delta', 'R22', 'R60', 'R124', 'R260', 'R1220', 'R1740']
# list repressors in order
repressors = [0, 0, 22, 60, 124, 260, 1220, 1740]
# add IPTG
IPTG_uM = [0, 0.1, 5, 10, 25, 50, 75, 100, 250, 500, 1000]

# initialize value with nan
df['strain'] = ['nan'] * len(df)
df['repressors'] = [0] * len(df)
df['IPTG_uM'] = [0] * len(df)

# Loop through rows in the plate to indicate the strain
for i, row in enumerate(list(string.ascii_uppercase[0:8])):
    # Extract index
    idx = [s for s in df.index if row in s]
    # indicate the corresponding strain
    df.set_value(idx, 'strain', strains[i]) 
    df.set_value(idx, 'repressors', repressors[i]) 
    df.set_value(idx, 'IPTG_uM', IPTG_uM) 

# Initialize pandas series to save the corrected YFP value
yfp_bgcorr = np.array([])

# Correct for the autofluorescence background
for i in np.arange(len(df)):
    data = df.iloc[i]
    auto = df[(df.IPTG_uM == data.IPTG_uM) &
              (df.strain == 'auto')].YFP_blank
    yfp_bgcorr = np.append(yfp_bgcorr, data.YFP_blank - auto)

yfp_bgcorr = pd.Series(yfp_bgcorr)
yfp_bgcorr.name = 'yfp_bgcorr'
yfp_bgcorr.index = df.index
df = pd.concat([df, yfp_bgcorr], join_axes=[df.index],
               axis=1, join='inner')

yfp_fc = np.array([])
# Compute the fold-change
for i in np.arange(len(df)):
    data = df.iloc[i]
    delta = df[(df.IPTG_uM == data.IPTG_uM) &
               (df.strain == 'delta')].yfp_bgcorr
    yfp_fc = np.append(yfp_fc, data.yfp_bgcorr / delta)

# Convert the fold-change to a pandas DataFrame.
yfp_fc = pd.Series(yfp_fc)
yfp_fc.name = 'fold_change'
yfp_fc.index = df.index
df = pd.concat([df, yfp_fc], join_axes=[df.index], axis=1, join='inner')

#============================================================================== 
# Plot all raw data
#============================================================================== 
rbs = df.strain.unique()
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
plt.savefig('output/' + str(date) + '_' + operator + \
            '_IPTG_titration_data.png')

#============================================================================== 
# Plot the WT control with the parameters determined from the MWC induction
# project
#============================================================================== 
# Define the IPTG concentrations to evaluate
IPTG = np.logspace(-7, -2, 100)
IPTG_lin = np.array([0, 1E-7])

energies = {'O1': -15.3, 'O2': -13.9, 'O3': -9.7, 'Oid': -17}
ea = -np.log(139.59)
ei = -np.log(0.53)


# Set the colors for the strains
colors = sns.color_palette('colorblind', n_colors=7)
colors[4] = sns.xkcd_palette(['dusty purple'])[0]

rep = np.sort(df[(df.strain!='auto') &\
        (df.strain!='delta')].repressors.unique())[::-1]
# Initialize figure
plt.figure()
for i, r in enumerate(rep):
    data = df[df.repressors==r]
    # Plot theoretical prediction
    # Log scale
    plt.plot(IPTG, mwc.fold_change_log(IPTG * 1E6,
        ea=ea, ei=ei, epsilon=4.5,
        R=r,
        epsilon_r=energies[data.operator.unique()[0]]),
        color=colors[i], label=None)
    # Linear scale
    plt.plot(IPTG_lin, mwc.fold_change_log(IPTG_lin * 1E6,
        ea=ea, ei=ei, epsilon=4.5,
        R=r,
        epsilon_r=energies[data.operator.unique()[0]]),
        linestyle='--', color=colors[i], label=None)

    # Plot data
    plt.scatter(data.IPTG_uM / 1E6, data.fold_change, color=colors[i], label=r)

plt.xscale('symlog', linthreshx=1E-7, linscalex=0.5)
plt.xlabel('IPTG (M)', fontsize=15)
plt.ylabel('fold-change', fontsize=16)
plt.ylim([-0.05, 1.1])
plt.xlim([-5E-9, np.max(IPTG)])
plt.legend(loc='upper left', title='repressors / cell')
plt.tick_params(labelsize=14)
plt.tight_layout()
plt.savefig('output/' + str(date) + '_' + operator + \
            '_theory_titration.png')

