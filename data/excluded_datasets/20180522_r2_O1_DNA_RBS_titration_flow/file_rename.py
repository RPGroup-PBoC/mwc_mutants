# -*- coding: utf-8 -*-
import numpy as np
import fcsparser
import os
import glob

# Define the details fo the expriment.
USERNAME = 'gchure'
DATE = 20180522
R = [0, 0, 60, 124, 260, 1220]
OPERATOR = 'O1'
FCS_PATTERN = 'RP2018-05-22'

savedir = '../../../data/flow/csv/'

# Define the order of rows and the cols.
ROWS = ('auto', 'delta', 'Y20I', 'Y20I', 'Y20I', 'Y20I')
COLS = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
RUN_NO = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
CLASS = ('DNA', 'DNA', 'DNA', 'DNA', 'DNA', 'DNA')

# Get the names of the files.
files = glob.glob('../../../data/flow/fcs/{}*.fcs'.format(FCS_PATTERN))
files = np.sort(files)

# Break the list up into columns.
ncols, nrows = len(COLS), len(ROWS)
col_groups = [files[i:i + nrows] for i in range(0, len(files), nrows)]
for i, col in enumerate(col_groups):
    for j, samp in enumerate(col):
        # Define the new name.
        name = '{0}_r{1}_{2}_R{3}_{4}_{5}uMIPTG'.format(
            DATE, RUN_NO[i], OPERATOR, R[j], ROWS[j], COLS[i])  

        # Load the file using fcsparser and save to csv.
        _, data = fcsparser.parse(samp)
        data.to_csv('{0}{1}.csv'.format(savedir, name))

        # Rename the fcs file.
        os.rename(samp, '../../../data/flow/fcs/{0}.fcs'.format(name))
