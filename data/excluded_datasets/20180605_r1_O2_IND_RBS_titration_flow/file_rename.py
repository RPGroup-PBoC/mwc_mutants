# -*- coding: utf-8 -*-
import numpy as np
import fcsparser
import os
import glob

# Define the details fo the expriment.
USERNAME = 'gchure'
DATE = 20180605
R = [0, 0, 260, 260, 260, 260, 260, 260]
OPERATOR = 'O2'
FCS_PATTERN = 'RP2018-06-05'

savedir = '../../../data/flow/csv/'

# Define the order of rows and the cols.
ROWS = ('auto', 'delta', 'Q294K', 'Q294K', 
        'Q294K', 'Q294K', 'Q294K','Q294K')
COLS = (0, 0.1, 5, 10, 25, 50, 75, 100, 250, 500, 1000, 5000)
RUN_NO = (1, 1, 1, 2, 3, 4, 5, 6)
CLASS = ('IND','IND','IND','IND','IND','IND')

# Get the names of the files
files = glob.glob('../../../data/flow/fcs/{0}*.fcs'.format(FCS_PATTERN))
files = np.sort(files)

# Break the list up into columns.
ncols, nrows = len(COLS), len(ROWS)
col_groups = [files[i:i + nrows] for i in range(0, len(files), nrows)]
for i, col in enumerate(col_groups):
    for j, samp in enumerate(col):
        # Define the new name.
        name = '{0}_r{1}_{2}_R{3}_{4}_{5}uMIPTG'.format(
            DATE, RUN_NO[j], OPERATOR, R[j], ROWS[j], COLS[i])  

        # Load the file using fcsparser and save to csv.
        _, data = fcsparser.parse(samp)
        data.to_csv('{0}{1}.csv'.format(savedir, name))

        # Rename the fcs file.
        os.rename(samp, '../../../data/flow/fcs/{0}.fcs'.format(name))
