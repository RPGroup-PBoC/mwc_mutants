import os
import glob
# Our numerical workhorses
import numpy as np
import pandas as pd
import scipy
# Import matplotlib stuff for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Seaborn, useful for graphics
import seaborn as sns

# favorite Seaborn settings for notebooks
rc={'lines.linewidth': 2, 
    'axes.labelsize' : 16, 
    'axes.titlesize' : 18,
    'axes.facecolor' : 'F4F3F6',
    'axes.edgecolor' : '000000',
    'axes.linewidth' : 1.2,
    'xtick.labelsize' : 13,
    'ytick.labelsize' : 13,
    'grid.linestyle' : ':',
    'grid.color' : 'a6a6a6'}
sns.set_context('notebook', rc=rc)
sns.set_style('darkgrid', rc=rc)
sns.set_palette("deep", color_codes=True)

# Import the project utils
import sys
sys.path.insert(0, '../../analysis/')

import mwc_induction_utils_processing as mwc
#=============================================================================== 
# define variables to use over the script
date = 20160928
username = 'sbarnes'
run = 'r1'

# define the patterns in the file names to read them
operator = 'O2'


# list the directory with the data
datadir = '../../../data/flow/csv/'
files = np.array(os.listdir(datadir))
csv_bool = np.array([str(date) in f and 'csv' in f for f in files])
files = files[np.array(csv_bool)]
mutants = np.array(['wt', 'Q294V', 'Q294K', 'F164T', 'Y20I', 'Q21M'])
rbs_wt = np.array(['auto', 'delta', 'RBS1027'])
rbs = 'RBS1027'
#repressors = np.array([0, 0, 870, 610, 130, 62, 30, 11])

concentrations = [0, 0.1, 5, 10, 25, 50, 75, 100, 250, 500, 1000, 5000] # uM IPTG

#=============================================================================== 
# define the parameter alpha for the automatic gating
alpha = 0.40

# initialize the DataFrame to save the mean expression levels
df = pd.DataFrame()
# read the files and compute the mean YFP value
for i, c in enumerate(concentrations):
    for j, strain in enumerate(mutants):
        # find the file
        try:
            if strain == 'wt':
                for rbs in rbs_wt:
                    r_file = glob.glob(datadir + str(date) + '_' + run + '_' + strain + '_' + \
                            operator + '_' + rbs + '_' + str(c) + 'uM' + '*csv')
                    print(r_file)
                    # read the csv file
                    dataframe = pd.read_csv(r_file[0])
                    # apply an automatic bivariate gaussian gate to the log front
                    # and side scattering
                    data = mwc.auto_gauss_gate(dataframe, alpha, 
                                                x_val='FSC-A', y_val='SSC-A',
                                                log=True)
                    # compute the mean and append it to the data frame along the
                    # operator and strain
                    df = df.append([[date, username, operator, strain, 
                                rbs, c,
                                data['FITC-A'].mean()]],
                                ignore_index=True)
            else:
                r_file = glob.glob(datadir + str(date) + '_' + run + '_' + strain + '_' + \
                        operator + '_' + rbs + '_' + str(c) + 'uM' + '*csv')
                print(r_file)
                # read the csv file
                dataframe = pd.read_csv(r_file[0])
                # apply an automatic bivariate gaussian gate to the log front
                # and side scattering
                data = mwc.auto_gauss_gate(dataframe, alpha, 
                                            x_val='FSC-A', y_val='SSC-A',
                                            log=True)
                # compute the mean and append it to the data frame along the
                # operator and strain
                df = df.append([[date, username, operator, strain, 
                            rbs, c,
                            data['FITC-A'].mean()]],
                            ignore_index=True)
        except:
            pass

# rename the columns of the data_frame
df.columns = ['date', 'username', 'operator', 'mutant', \
        'rbs', 'IPTG_uM', 'mean_YFP_A']

# initialize pandas series to save the corrected YFP value
mean_bgcorr_A = np.array([])
# correct for the autofluorescence background
for i in np.arange(len(df)):
    data = df.loc[i]
    auto = df[(df.IPTG_uM == data.IPTG_uM) & \
              (df.rbs == 'auto')].mean_YFP_A
    mean_bgcorr_A = np.append(mean_bgcorr_A, data.mean_YFP_A - auto)

mean_bgcorr_A = pd.Series(mean_bgcorr_A)
mean_bgcorr_A.name = 'mean_YFP_bgcorr_A'
df = pd.concat([df, mean_bgcorr_A], join_axes=[df.index],
                axis=1, join='inner')

mean_fc_A = np.array([])
# compute the fold-change
for i in np.arange(len(df)):
    data = df.loc[i]
    delta = df[(df.IPTG_uM == data.IPTG_uM) & \
              (df.rbs == 'delta')].mean_YFP_bgcorr_A
    mean_fc_A = np.append(mean_fc_A, data.mean_YFP_bgcorr_A / delta)

mean_fc_A = pd.Series(mean_fc_A)
mean_fc_A.name = 'fold_change_A'
df = pd.concat([df, mean_fc_A], join_axes=[df.index], axis=1, join='inner')

# write
df.to_csv('output/' + str(date) + '_' + run + '_' + operator + \
        '_IPTG_titration_MACSQuant.csv', index=False)
#=============================================================================== 
# Add the comments to the header of the data file
filenames = ['./comments.txt', 'output/' + str(date) + '_' + run + '_' + \
             operator + '_IPTG_titration_MACSQuant.csv']
with open('../../../data/' + str(date) + '_' + run + '_' + operator + \
        '_IPTG_titration_MACSQuant.csv', 'w') as output:
    for fname in filenames:
        with open(fname) as infile:
            output.write(infile.read())
