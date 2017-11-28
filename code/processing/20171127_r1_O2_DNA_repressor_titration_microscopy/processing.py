import numpy as np
import matplotlib.pyplot as plt
import pboc.plotting
import pboc.image
import skimage.io
import pandas as pd
import glob
import seaborn as sns


# Define the constants of the experiment.
DATE = 20171127
OPERATOR = 'O2'
RUN_NO = 1
USERNAME = 'mrazomej'
MUT_CLASS = 'DNA'
IP_DIST = 0.16  # in units of µm per pixel

# Load the images of the slides.
slide_files = glob.glob('../../../data/images/*profile*/*/*.tif')
noise_files = glob.glob('../../../data/images/*noise*/*/*.tif')
slide_ims = skimage.io.ImageCollection(slide_files)
noise_ims = skimage.io.ImageCollection(noise_files)
mean_field = pboc.image.projection(slide_ims, mode='mean')
mean_dark = pboc.image.projection(noise_ims, mode='mean')

# Load the images.
files = glob.glob('../../../data/images/{0}*IPTG*/*/'.format(DATE))

# Choose a random file to generate an example segmentation.
rand_pos = np.random.choice(files)
columns = ['date', 'username', 'mut_class', 'mutant', 'rbs', 'repressors',
           'mean_intensity', 'area', 'eccentricity', 'operator', 'IPTG_uM']

df = pd.DataFrame([], columns=columns)
for i, f in enumerate(files):
    f
    # Get the identifying information
    date, mut, op, rbs, IPTG, _ = f.split('/')[-3].split('_')
    IPTG = float(IPTG.split('uM')[0])
    if (rbs != 'auto') and (rbs != 'delta'):
        rep = int(rbs.split('R')[-1])
    else:
        rep = 0

    # Load the three images.
    phase = skimage.io.imread('{0}img_000000000_Brightfield_000.tif'.format(f))
    vol = skimage.io.imread('{0}img_000000000_TRITC_000.tif'.format(f))
    yfp = skimage.io.imread('{0}img_000000000_YFP_000.tif'.format(f))

    # Segment the volume marker image.
    seg = pboc.image.log_segmentation(vol, label=True)

    # Correct the intensity image.
    yfp_flat = pboc.image.generate_flatfield(yfp, mean_dark, mean_field)

    # Measure the properties and store them.
    props = skimage.measure.regionprops(seg, yfp_flat)
    for prop in props:
        area = prop.area * IP_DIST**2
        eccentricity = prop.eccentricity
        mean_intensity = prop.mean_intensity
        cell_dict = dict(date=DATE, username=USERNAME, mut_class=MUT_CLASS,
                         mutant=mut, repressors=rep,
                         mean_intensity=mean_intensity, area=area,
                         eccentricity=eccentricity, operator=OPERATOR,
                         IPTG_uM=IPTG, rbs=rbs)
        df = df.append(cell_dict, ignore_index=True)

    # Make the example segmentation if specified.
    if f == rand_pos:
        BAR_LENGTH = 10
        # Make a merge.
        phase_float = (phase - phase.min()) / (phase.max() - phase.min())
        phase_float[20:30, 20:20 + int(BAR_LENGTH / IP_DIST)] = 1.0
        phase_copy = np.copy(phase_float)
        phase_copy[seg > 0] = 0.8
        merge = np.dstack((phase_float, phase_copy, phase_copy))
        savename = '{0}/{1}_{2}_example_segmentation.png'.format('output', DATE,
                                                                 MUT_CLASS)
        with sns.axes_style('white'):
            fig, ax = plt.subplots(1, 1)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_frame_on(False)
            ax.imshow(merge)
            plt.tight_layout()
            plt.savefig(savename, bbox_inches='tight')


#%% Prune the data frame.
area_bounds = [0.5, 6.0]
ecc_bound = 0.8

data = df[(df['area'] >= area_bounds[0]) & (df['area'] <= area_bounds[1])
          & (df['eccentricity'] > ecc_bound)]

# Save it to disk.
target = '{0}/{1}_{2}_{3}_repressor_titration_microscopy_measurements.csv'.format('output',
                                                                                  DATE, OPERATOR,
                                                                                  MUT_CLASS)

data.to_csv(target, index=False)
