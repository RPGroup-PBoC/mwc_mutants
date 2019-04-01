# -*- coding: utf-8 -*-
import numpy as np
import fcsparser
import os
import glob

# Define the details fo the expriment.
USERNAME = 'gchure'
DATE = 20190326
RUN_NO = 1 
FCS_PATTERN = 'RP2019-03-26'

savedir = '../../../data/flow/csv/'

# Define the order of rows and the cols.
R = (0, 0, 260, 260, 260, 260)
ROWS = ('auto', 'delta', 'F164T', 'Q294V', 'Q294K', 'Q294R')
OPS = ('NA', 'O3', 'O3', 'O3', 'O3', 'O3')
COLS = (0, 0.1, 5, 10, 25, 50, 75, 100, 250, 500, 1000, 5000)

# Get the names of the files
files = glob.glob('../../../data/flow/fcs/{0}*r{1}*.fcs'.format(FCS_PATTERN, RUN_NO))
files = np.sort(files)

# Break the list up into columns.
ncols, nrows = len(COLS), len(ROWS)
col_groups = [files[i:i + nrows] for i in range(0, len(files), nrows)]
for i, col in enumerate(col_groups):
    for j, samp in enumerate(col):
        # Define the new name.
        name = '{0}_r{1}_{2}_R{3}_{4}_{5}uMIPTG'.format(
            DATE, RUN_NO, OPS[j], R[j], ROWS[j], COLS[i])  

        # Load the file using fcsparser and save to csv.
        _, data = fcsparser.parse(samp)
        data.to_csv('{0}{1}.csv'.format(savedir, name))

        # Rename the fcs file.
        os.rename(samp, '../../../data/flow/fcs/{0}.fcs'.format(name))