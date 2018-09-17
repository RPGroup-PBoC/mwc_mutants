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
# Plot the theory vs data for all 4 operators with the credible region
#===============================================================================
# Define the IPTG concentrations to evaluate
IPTG = np.logspace(-7, -2, 100)
IPTG_lin = np.array([0, 1E-7])

# Group by class of mutant, operator and repressor copy number
df_group = df.groupby(['class', 'operator', 'repressors'])

# Extract the available groups
groups = [g[0] for g in df_group]

# Save plots in multipage PDF
with PdfPages(output + 'theory_vs_data_epsilonRA.pdf') as pdf:
    for group, data in df_group:
        # Define color palette
        colors = mwc.color_dict(np.append(data.strain.unique(), 'wt'), 'colorblind')
        # Initialize the plot to set the size
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        # group data by unique strains
        # extract the wt data for those sets
        dates = data.date.unique()
        df_wt = df[(df.date.isin(dates)) & (df.strain=='wt')\
                & (df.operator==data.operator.unique()[0])]
        # Load pre-computed fold-change with credible region
        df_fc = pd.read_csv(datadir + 'mcmc/epsilonRA/fold_change/' +\
                'fold_change_' + df_wt.operator.unique()[0] +\
                '_' + str(df_wt.repressors.unique()[0]) + '_wt.csv')
        # sort by IPTG concentration
        df_fc.sort_values('IPTG_M', inplace=True)
        # plot the theory
        # Log-scale
        ax.plot(df_fc.IPTG_M[df_fc.IPTG_M >= 1E-7], 
                df_fc.fold_change_map[df_fc.IPTG_M >= 1E-7],
                zorder=1, color=colors['wt'], label='')
        # Linear scale
        ax.plot(df_fc.IPTG_M[df_fc.IPTG_M <= 1E-7], 
                df_fc.fold_change_map[df_fc.IPTG_M <= 1E-7],
                zorder=1, linestyle='--', color=colors['wt'], label='')
        # Plot the credible region
        ax.fill_between(df_fc.IPTG_M, df_fc.hpd_min, df_fc.hpd_max,
                        alpha=0.3, color=colors['wt'])

        # compute the mean value for each concentration
        fc_mean = df_wt.groupby('IPTG_uM').fold_change.mean()
        # compute the standard error of the mean
        fc_err = df_wt.groupby('IPTG_uM').fold_change.std() / \
        np.sqrt(df_wt.groupby('IPTG_uM').size())
        ax.errorbar(np.sort(df_wt.IPTG_uM.unique()) / 1E6,
                fc_mean, yerr=fc_err, fmt='o', 
                label= 'wild-type', zorder=100, color=colors['wt'])

        # Plot the mutants data
        data_strain = data.groupby('strain')
        for strain, sdata in data_strain:
            # Load pre-computed fold-change with credible region
            df_fc = pd.read_csv(datadir + 'mcmc/epsilonRA/fold_change/' +\
                    'fold_change_' + sdata.operator.unique()[0] +\
                    '_' + str(sdata.repressors.unique()[0]) + '_' + strain + '.csv')
            # sort by IPTG concentration
            df_fc.sort_values('IPTG_M', inplace=True)
            # plot the theory
            # Log-scale
            ax.plot(df_fc.IPTG_M[df_fc.IPTG_M >= 1E-7], 
                    df_fc.fold_change_map[df_fc.IPTG_M >= 1E-7],
                    zorder=1, color=colors[strain], label='')
            # Linear scale
            ax.plot(df_fc.IPTG_M[df_fc.IPTG_M <= 1E-7], 
                    df_fc.fold_change_map[df_fc.IPTG_M <= 1E-7],
                    zorder=1, linestyle='--', color=colors[strain], label='')
            ax.fill_between(df_fc.IPTG_M, df_fc.hpd_min, df_fc.hpd_max,
                        alpha=0.3, color=colors[strain])

            # compute the mean value for each concentration
            fc_mean = sdata.groupby('IPTG_uM').fold_change.mean()
            # compute the standard error of the mean
            fc_err = sdata.groupby('IPTG_uM').fold_change.std() / \
            np.sqrt(sdata.groupby('IPTG_uM').size())
            ax.errorbar(np.sort(sdata.IPTG_uM.unique()) / 1E6,
                    fc_mean, yerr=fc_err, fmt='o', 
                    label=strain, zorder=100, color=colors[strain])
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
