import os
import glob
import shutil
import numpy as np


# Define the experimental constants.
DATE = 20171122
OPERATOR = 'O2'

CORRECTION_MUTS = {'Y20IREAL': 'Y20I', 'Y20I': 'Q21A'}
CORRECTION_R = 'R1220'


# Get the names of all of the files.
files = glob.glob(
    '../../../../data/images/incorrectly_named/{0}*IPTG*'.format(DATE))

files
for i, f in enumerate(files):
    date, mut, op, rep, IPTG, num = f.split('/')[-1].split('_')
    if (mut in CORRECTION_MUTS.keys()) and (rep == CORRECTION_R):
        mut = CORRECTION_MUTS[mut]

    # Define the new name.
    new_name = '{0}_{1}_{2}_{3}_{4}_{5}'.format(date, mut, op, rep, IPTG, num)

    # Move it to the new place.
    shutil.copytree(f, '{0}/{1}'.format('../../../../data/images', new_name))

print('I undid your fuck up. Thank you -- come again')
