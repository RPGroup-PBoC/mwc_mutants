# Import dependencies.
import os
import sys
import itertools
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
import sys
sys.path.insert(0, '../../')

import mwc_mutants_utils as mwc
mwc.set_plotting_style()

#============================================================================== 
# Variables to edit on each script
#============================================================================== 
# Define variables to use over the script
date = 20170308
username = 'mrazomej'
run = 'r2'
operator = 'O2'
mutant_class = 'dna' # options dna/inducer/double
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

    # List all strains from the project
    strain_no_mut = ['wt', 'delta', 'auto']
    # DNA binding site mutants
    strain_dna = ['Y20I', 'Q21M', 'Q21A']
    # Inducer binding site
    strain_inducer = ['F164T', 'Q294K', 'Q294V']
    # Double mutants
    strain_double = [d[0] + '-' + d[1] for d in itertools.product(strain_dna,
                     strain_inducer)] 
    # Extract strain
    strain = list(set(f_split).intersection(strain_no_mut + strain_dna +\
                                            strain_inducer))[0]

    # Infer the class of mutant. This can be nan, dna, inducer or double
    if strain in strain_no_mut:
        mut_class = 'NA'
    elif strain in strain_dna:
        mut_class = 'dna'
    elif strain in strain_inducer:
        mut_class = 'inducer'
    elif strain in strain_double:
        mut_class = 'double'

    return dict(zip(['op', 'rep', 'iptg', 'strain', 'class'],
                    [op, rep, iptg, strain, mut_class]))

#============================================================================== 

# List the target directory.
datadir = '../../../data/flow/csv/'
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
    print(filename)
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
                    file_info['op'], file_info['class'],
                    file_info['strain'], file_info['rep'],
                    file_info['iptg'], data['FITC-A'].mean()]],
                    ignore_index=True)

# Rename the columns of the data_frame
df.columns = ['date', 'username', 'operator', 'class', 'strain',
              'repressors', 'IPTG_uM', 'mean_YFP']

# Initialize pandas series to save the corrected YFP value
mean_bgcorr = np.array([])

# Correct for the autofluorescence background
for i in np.arange(len(df)):
    data = df.loc[i]
    auto = df[(df.IPTG_uM == data.IPTG_uM) &
              (df.strain == 'auto')].mean_YFP
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
               (df.strain == 'delta')].mean_YFP_bgcorr
    mean_fc = np.append(mean_fc, data.mean_YFP_bgcorr / delta)

# Convert the fold-change to a pandas DataFrame.
mean_fc = pd.Series(mean_fc)
mean_fc.name = 'fold_change'
df = pd.concat([df, mean_fc], join_axes=[df.index], axis=1, join='inner')

# write
df.to_csv('output/tmp.csv', index=False)

# Add the comments to the header of the data file
filenames = ['./comments.txt', 'output/tmp.csv']

with open('output/' + str(date) + '_' + run + '_' + operator + '_' + \
          mutant_class + '_IPTG_titration_MACSQuant.csv', 'w') as output:
    for fname in filenames:
        with open(fname) as infile:
            output.write(infile.read())

# Remove temporary file
os.remove(filenames[1])
