# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mut.thermo
import mut.stats
import mut.viz
color = mut.viz.color_selector('mut')
pboc = mut.viz.color_selector('pboc')
constants = mut.thermo.load_constants()
mut.viz.set_plotting_style()

# Load the double data
data = pd.read_csv('../../data/csv/summarized_data.csv')
predictions = pd.read_csv('../../data/csv/Fig4_O2_double_samples.csv')

predictions
# Load the necessary statistics samples. 
