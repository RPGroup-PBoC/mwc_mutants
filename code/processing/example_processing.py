"""
Title:
    example_processing.py
Author:
    Griffin Chure and Manuel Razo-Mejia
Creation Date:
    20170220
Last Modified:
    20170220
Purpose:
    This script serves as a representative example of our data processing
    pipeline. This script reads in a set of csv files containing the output
    from the MACSQuant Flow Cytomter (after being converted from the flow
    cytometry standard .fcs format), peforms unsupervised gating, and computes
    the measured fold-change in gene expression.
"""

# Import dependencies.
import os
import glob
import re
import numpy as np
import pandas as pd
import scipy

# Import matplotlib stuff for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Seaborn, useful for graphics
import seaborn as sns

# Set the plotting style.
import mwc_induction_utils as mwc
mwc.set_plotting_style()

#============================================================================== 
# Variables to edit on each script
#============================================================================== 
# Define variables to use over the script
date = 20170101
username = 'username'
run = 'rx'
operator = 'Ox'
mutant = 'dna' # options are dna/inducer/double
#============================================================================== 
# Define function that extracts information from file name
def file_info_parse(file):
    '''
    Extracts the relevant information from the file name by splitting the
    string
    '''
    # Split the file name by the underscore
    f_split = re.split(pattern='_', string=file)

    # Find the operator
    op_regex = re.compile('O.')
    op = [f for f in f_split if re.match(op_regex, f)][0]

    # Find the repressor copy number
    rep_regex = re.compile('R.')
    rep = [f for f in f_split if re.match(rep_regex, f)][0]
    rep = int(re.findall('\d+', rep)[0])

    # Find the IPTG concentration
    iptg_str = [f for f in f_split if 'uM' in f][0]
    iptg = re.findall('\d+\.\d+', iptg_str)
    if not iptg:
        iptg = re.findall('\d+', iptg_str)
    iptg = float(iptg[0])

    # Extract mutant or auto/delta
    mut = f_split[2]

    return dict(zip(['operator', 'repressor', 'IPTG_uM', 'mutant'],
                    [op, rep, iptg, mut]))

#============================================================================== 

# List the target directory.
datadir = '../../data/'
files = np.array(os.listdir(datadir))
# Find files that contain the date, run number and are CSV
csv_bool = np.array([str(date) in f and run in f \
                     and 'csv' in f for f in files])
files = files[np.array(csv_bool)]


# Define the parameter alpha for the automatic gating
alpha = 0.40

# Initialize the DataFrame to save the mean expression levels
df = pd.DataFrame()
# Read the files and compute the mean YFP value
for filename in files:
    file_info = file_info_parse(filename)
    dataframe = pd.read_csv(datadir + filename)
    # Apply an automatic bivariate gaussian gate to the log front
    # and side scattering
    data = mwc.auto_gauss_gate(dataframe, alpha,
                                x_val='FSC-A', y_val='SSC-A',
                                log=True)
    # Compute the mean and append it to the data frame along the
    # operator and strain
    df = df.append([[date, username, 
                    file_info['op'], file_info['mut'], 
                    file_info['rep'], file_info['iptg'],
                    data['FITC-A'].mean()]],
                    ignore_index=True)

# Rename the columns of the data_frame
df.columns = ['date', 'username', 'operator', 'repressors', 'IPTG_uM', 
              'mean_YFP']

# Initialize pandas series to save the corrected YFP value
mean_bgcorr = np.array([])

# Correct for the autofluorescence background
for i in np.arange(len(df)):
    data = df.loc[i]
    auto = df[(df.IPTG_uM == data.IPTG_uM) &
              (df.rbs == 'auto')].mean_YFP
    mean_bgcorr = np.append(mean_bgcorr, data.mean_YFP - auto)

mean_bgcorr = pd.Series(mean_bgcorr)
mean_bgcorr.name = 'mean_YFP_bgcorr'
df = pd.concat([df, mean_bgcorr], join_axes=[df.index],
               axis=1, join='inner')
mean_fc = np.array([])

# Compute the fold-change
for i in np.arange(len(df)):
    data = df.loc[i]
    delta = df[(df.IPTG_uM == data.IPTG_uM) &
               (df.rbs == 'delta')].mean_YFP_bgcorr
    mean_fc = np.append(mean_fc, data.mean_YFP_bgcorr / delta)

# Convert the fold-change to a pandas DataFrame.
mean_fc = pd.Series(mean_fc)
mean_fc.name = 'fold_change'
df = pd.concat([df, mean_fc], join_axes=[df.index], axis=1, join='inner')

# write
df.to_csv('output/' + str(date) + '_' + run + '_' + operator + \
        mutant + '_IPTG_titration_MACSQuant.csv', index=False)

# Add the comments to the header of the data file
filenames = ['./comments.txt', 'output/tmp.csv']

with open('output/' + str(date) + '_' + run + '_' + operator + \
          mutant + '_IPTG_titration_MACSQuant.csv', 'w') as output:
    for fname in filenames:
        with open(fname) as infile:
            output.write(infile.read())

# Remove temporary file
os.remove(filenames[1])
