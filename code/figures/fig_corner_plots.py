import os
import glob
import pickle
import re
import glob

# Our numerical workhorses
import numpy as np
import pandas as pd

# Import the project utils
import sys
sys.path.insert(0, '../')
import mwc_mutants_utils as mwc

# Import matplotlib stuff for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages

# Seaborn, useful for graphics
import seaborn as sns

# Library for corner plots
import corner

mwc.set_plotting_style()

#===============================================================================
# Set output directory based on the graphicspath.tex file to print in dropbox
#===============================================================================
dropbox = open('dropbox_path.txt')
output = dropbox.read()
output = output[:-1]

#===============================================================================
# List the MCMC flat chains
#===============================================================================
# read the flow_master
datadir = '../../data/mcmc/'

files = glob.glob(datadir + '*pkl')

#===============================================================================
# Plot available flat chains
#===============================================================================
# Save plots in multipage PDF
with PdfPages(output + 'corner_plots.pdf') as pdf:
    for i, filename in enumerate(files):
        print(filename)
        # Initialize the plot to set the size
        with open('../../data/mcmc/' + filename, 'rb') as file:
            unpickler = pickle.Unpickler(file)
            gauss_flatchain = unpickler.load()
            gauss_flatlnprobability = unpickler.load()

        # Draw the corner plot
        chain = np.vstack([gauss_flatchain[:,0], np.exp(-gauss_flatchain[:,1]),
                        np.exp(-gauss_flatchain[:,2])]).T
        fig, ax = plt.subplots(3, 3, figsize=(7, 7))
        corner.corner(chain, bins=50, plot_contours=True,
                        labels=[r'$\Delta\varepsilon_{RA}$', r'$K_A$', 
                                r'$K_I$'], fig=fig)
        plt.suptitle(re.split(pattern='/', string=filename)[-1])
        plt.tight_layout()
        pdf.savefig()
        plt.close()
