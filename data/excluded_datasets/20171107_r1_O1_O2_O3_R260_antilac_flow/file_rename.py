import numpy as np
import fcsparser
import os
import glob

# Define the details fo the expriment.
USERNAME = 'gchure'
DATE = 20171107
MUTANT = 'antilac'
CLASS = 'other'
R = 260
RUN_NO = 1
FCS_PATTERN = 'RP2017'

savedir = '../../../data/flow/csv/'
# Define the order of rows and the cols.
ROWS = ('auto', 'O1_delta', 'O2_delta', 'O3_delta', 'O1_RBS1027', 'O2_RBS1027',
        'O3_RBS1027')

COLS = (0, 0.1, 5, 10, 25, 50, 75, 100, 250, 500, 1000, 5000)


# Get the names of the files.
files = glob.glob('../../../data/flow/fcs/{0}*.fcs'.format(FCS_PATTERN))
files = np.sort(files)

# Break the list up into columns.
ncols, nrows = len(COLS), len(ROWS)
col_groups = [files[i:i + nrows] for i in range(0, len(files), nrows)]
for i, col in enumerate(col_groups):
    for j, samp in enumerate(col):
        # Define the new name.
        name = '{0}_r{1}_{2}_{3}_{4}uMIPTG'.format(
            DATE, RUN_NO, MUTANT, ROWS[j], COLS[i])

        # Load the file using fcsparser and save to csv.
        _, data = fcsparser.parse(samp)
        data.to_csv('{0}{1}.csv'.format(savedir, name))

        # Rename the fcs file.
        os.rename(src=samp, dst='../../../data/flow/fcs/{0}.fcs'.format(name))
