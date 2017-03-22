import os
import glob
import pickle
import re

# Our numerical workhorses
import numpy as np
import pandas as pd

# Import the project utils
import sys
sys.path.insert(0, '../')
import mwc_mutants_utils as mwc

# Import matplotlib stuff for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages

# Seaborn, useful for graphics
import seaborn as sns

mwc.set_plotting_style()

#===============================================================================
# Set output directory based on the graphicspath.tex file to print in dropbox
#===============================================================================
dropbox = open('dropbox_path.txt')
output = dropbox.read()
output = output[:-1]

#===============================================================================
# Read the data
#===============================================================================
# read the flow_master
datadir = '../../data/'
df = pd.read_csv(datadir + 'flow_master.csv', comment='#', index_col=0)

# Now we remove the autofluorescence and delta values
df = df[(df.strain != 'auto') & (df.strain != 'delta')]

#===============================================================================
# Plot available raw data
#===============================================================================
# Define the IPTG concentrations to evaluate
IPTG = np.logspace(-7, -2, 100)
IPTG_lin = np.array([0, 1E-7])

energies = {'O1': -15.3, 'O2': -13.9, 'O3': -9.7, 'Oid': -17}
ka = -np.log(139.59)
ki = -np.log(0.53)

# Set the colors for the strains
colors = sns.color_palette('colorblind', n_colors=7)
colors[4] = sns.xkcd_palette(['dusty purple'])[0]

# Group by class of mutant, operator and repressor copy number
df_group = df.groupby(['class', 'operator', 'repressors'])

# Extract the available groups
groups = [g[0] for g in df_group]

# Save plots in multipage PDF
with PdfPages(output + 'raw_data.pdf') as pdf:
    for group, data in df_group:
        # Initialize the plot to set the size
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        # group data by unique strains
        # extract the wt data for those sets
        dates = data.date.unique()
        df_wt = df[(df.date.isin(dates)) & (df.strain=='wt')]
        # compute the mean value for each concentration
        fc_mean = df_wt.groupby('IPTG_uM').fold_change.mean()
        # compute the standard error of the mean
        fc_err = df_wt.groupby('IPTG_uM').fold_change.std() / \
        np.sqrt(df_wt.groupby('IPTG_uM').size())
        ax.errorbar(np.sort(df_wt.IPTG_uM.unique()) / 1E6,
                fc_mean, yerr=fc_err, fmt='o', 
                label= 'wild-type', zorder=100)
        # Plot theoretical prediction
        # Log scale
        ax.plot(IPTG, mwc.fold_change_log(IPTG * 1E6,
            ka=ka, ki=ki, epsilon=4.5,
            R=df_wt.repressors.unique()[0],
            epsilon_RA=energies[df_wt.operator.unique()[0]]),
            color='black', label='prediction')
        # Linear scale
        ax.plot(IPTG_lin, mwc.fold_change_log(IPTG_lin * 1E6,
            ka=ka, ki=ki, epsilon=4.5,
            R=df_wt.repressors.unique()[0],
            epsilon_RA=energies[df_wt.operator.unique()[0]]),
            linestyle='--', color='black', label=None)

        # Plot the mutants data
        data_strain = data.groupby('strain')
        for strain, sdata in data_strain:
            # compute the mean value for each concentration
            fc_mean = sdata.groupby('IPTG_uM').fold_change.mean()
            # compute the standard error of the mean
            fc_err = sdata.groupby('IPTG_uM').fold_change.std() / \
            np.sqrt(sdata.groupby('IPTG_uM').size())
            ax.errorbar(np.sort(sdata.IPTG_uM.unique()) / 1E6,
                    fc_mean, yerr=fc_err, fmt='o', 
                    label= strain, zorder=100)
            ax.set_xscale('symlog', linthreshx=1E-7, linscalex=0.5)
            ax.set_xlabel('IPTG (M)', fontsize=15)
            ax.set_ylabel('fold-change', fontsize=16)
            ax.set_ylim([-0.01, 1.1])
            ax.set_xlim(left=-5E-9)
            ax.tick_params(labelsize=14)
            ax.legend(loc='upper left')
            ax.set_title(str(group))
        plt.tight_layout()
        pdf.savefig()
        plt.close()
