#%%
# -*- coding: utf-8 -*-
import numpy as np
import sys
import glob
import matplotlib.pyplot as plt
import pandas as pd
import skimage.io
import skimage.segmentation
import skimage.morphology
sys.path.insert(0, '../../../')
import mut.image
import mut.viz
import mut.stats
colors = mut.viz.pub_style()


# Define the experimental parameters
DATE = 20180404
USERNAME = 'gchure'
RUN_NO = 2
CLASS = 'DNA'
MUT = ['Q21M']
OP = 'O2'
IPTG = [0]
IP_DIST = 0.16 # in µm per pixel

# %%
data_dir = '../../../data/images/20180404_run2/'
# data_dir = 'data/images/'

# Load the images for the flat field
noise_files = glob.glob('{}/{}_camera_noise*/*.tif'.format(data_dir, DATE))
noise_ims = [skimage.io.ImageCollection(f)[0] for f in noise_files]
avg_noise = mut.image.projection(noise_ims, mode='mean')
yfp_files = glob.glob('{}/{}_YFP_profile_*/*.tif'.format(data_dir, DATE))
yfp_ims = [skimage.io.ImageCollection(f)[0] for f in yfp_files]
avg_field = mut.image.projection(yfp_ims, mode='mean')

# %% Segment auto fluorescence and store in a dataframe.
df = pd.DataFrame([], columns=['date', 'username', 'mutant', 'operator', 'strain', 
                               'IPTGuM', 'repressors', 'mean_intensity',
                               'area', 'eccentricity'])
samples = glob.glob('{}/{}_*/*.tif'.format(data_dir, DATE))
for i, s in enumerate(samples):
    if ('camera_noise' not in s) and ('YFP_profile') not in s:
        try: 
            _, mutant, strain, op, IPTG, _ = s.split('/')[-2].split('_')
            IPTG = float(IPTG.split('uM')[0])
            R = int(strain[1:])
            op = 'O2'
        except:
            _, strain, _ = s.split('/')[-2].split('_')
            R = 0
            IPTG = 0
            op = 'O2'
            mutant = 'none'

        # Median filter the fluorescence image.
        im = skimage.io.ImageCollection(s)

        seg = mut.image.log_segmentation(im[1], thresh=0.0001)
        seg_label = skimage.measure.label(seg)
        selem = skimage.morphology.square(3)

        im_filt = skimage.filters.median(im[2], selem)

        # Flatten the image. 
        im_flat = mut.image.generate_flatfield(im_filt, avg_noise, avg_field, median_filt=False )

        # Compute the properties and assemble the dataframe.
        props = skimage.measure.regionprops(seg_label, im_flat)
        for p in props:
            cell_dict = {'strain':strain, 'area':p.area * IP_DIST**2, 
                       'eccentricity':p.eccentricity,
                       'mean_intensity':p.mean_intensity,
                       'date': DATE, 'username':USERNAME, 'mutant':mutant,
                       'operator': op, 'IPTGuM':IPTG, 'repressors': R}
            df = df.append(cell_dict, ignore_index=True)

# %% Prune the data frame to restrict to single cells.
filt_df = df[(df['area'] >= 1) & (df['area'] <= 4) & (df['eccentricity'] >= 0.8)]

_ = plt.hist(df['eccentricity'], bins=100)
#%%  Compute the fold-change.
mean_auto = filt_df[filt_df['strain']=='auto']['mean_intensity'].mean()
mean_delta = filt_df[filt_df['strain']=='delta']['mean_intensity'].mean()
def compute_fc(df):
    fc = (df['mean_intensity'].mean() - mean_auto) / (mean_delta - mean_auto)
    new_df = {'date': df.date.unique()[0], 'mutant':df.mutant.unique()[0],
    'repressors':df.repressors.unique()[0], 'username':df.username.unique()[0], 'mean_YFP':df.mean_intensity.mean(),
    'fold_change': fc, 'run_number': RUN_NO}
    fc = pd.Series(new_df)
    return fc
grouped = filt_df.groupby(['strain']).apply(compute_fc)
fc_df = pd.DataFrame(grouped).reset_index()
fc_df.to_csv('./output/{}_r{}_{}_{}_microscopy_fold_change.csv'.format(DATE, RUN_NO, MUT[0], OP),
index=False)


# Plot the fold-change vs repressor copy number
fig, ax = plt.subplots(1, 1)
ax.set_xlabel('repressor copy number')
ax.set_ylabel('fold-change')
fc_df = fc_df[fc_df['repressors'] != 0]
fc_df.sort_values('repressors', inplace=True)
_ = ax.loglog(fc_df['repressors'], fc_df['fold_change'], '--o')
plt.tight_layout()
plt.savefig('./output/{}_r{}_{}_{}_microscopy_fold_change.png'.format(DATE, RUN_NO, MUT[0],
OP), bbox_inches='tight')


print(fc_df)