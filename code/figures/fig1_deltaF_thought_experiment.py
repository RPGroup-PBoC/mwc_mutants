# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
constants = mut.thermo
colors = mut.viz.color_selector('pboc')
mut.viz.plotting_style()



fig, ax = plt.subplots(1, 1, figsize=(2.5, 2.5))
ax.hlines(0, 1E-2, 1E4, 'k', linestyle=':')