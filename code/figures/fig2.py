# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.viz
pboc = mut.viz.color_selector('pboc')
color = mut.viz.color_selector('mut')
mut.viz.plotting_style()

# Load the summarized data. 
data = pd.read_csv('../../data/csv/summarized_data.csv')

# Instantiate the figure
fig, ax = plt.subplots(3, 2, figsize=(7, 6))

# Format the axes as necessary
ax[0,0].axis('off')