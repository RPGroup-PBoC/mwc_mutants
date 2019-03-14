import sys
sys.path.insert(0, '../../')
import mut.viz
import mut.thermo
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
import itertools
from bokeh.models import ColumnDataSource, Div
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select, Slider, RadioButtonGroup, Button
import glob
constants = mut.thermo.load_constants()
pboc = mut.viz.color_selector('pboc')

def induction_plot():
    data = pd.read_csv('../../data/csv/RazoMejia2018_data.csv')
    data.loc['repressors'] *= 2
    data.rename(columns={'fold_change_A':'fold_change',
                        'IPTG_uM':'IPTGuM'})
