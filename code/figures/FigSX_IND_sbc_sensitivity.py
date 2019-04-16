# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.viz
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()

# Load the sbc data
sbc_data = pd.read_csv('../../data/csv/IND_sbc.csv')



# ##############################################################################
# FIGURE INSTANTIATION 
# ##############################################################################
fig, ax = plt.subplots(2, 2)
