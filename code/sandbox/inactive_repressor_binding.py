import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import mut.thermo
import mut.viz
mut.viz.plotting_style()
colors = mut.viz.color_selector('pboc')


# Load the brewster and hernan data
data = pd.read_csv('../../data/csv/Garcia2011_Brewster2014_data.csv')

