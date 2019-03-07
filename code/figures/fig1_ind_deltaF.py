# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mut.thermo
import mut.viz
pboc = mut.viz.color_selector('pboc')

constants = mut.thermo.load_constants()
data = pd.read_csv('../../data/csv/RazoMejia2018_delta_F_statistics.csv')

rep_colors = {22:pboc['red'], 60:pboc['blue'], 124:pboc['green'],
             260:pboc['purple'], 1220:pboc['dark_green'], 1740:pboc['dark_brown']}

operator_glyphs = {'O1': 's', 'O2': 'o', 'O3': '^'}

# Define the reference states.
c0 = 50; # in µM
R0 = constants['RBS1027']
epRA0 = constants['O2']


