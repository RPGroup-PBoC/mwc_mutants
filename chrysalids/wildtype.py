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
import os
constants = mut.thermo.load_constants()
pboc = mut.viz.color_selector('pboc')

def induction_app(data):
    # #######################
    # INSTANTIATING FIGURES 
    # ####################### 
    p_O1 = bokeh.plotting.figure(width=400, height=300, 
                                x_axis_type='log',
                                y_axis_label='fold-change',
                                x_axis_label='IPTG [µM]',
                                title='    O1 | ∆εRA = -15.3 kT    ',
                                y_range=[-0.1, 1.1])
    p_O2 = bokeh.plotting.figure(width=400, height=300, 
                                x_axis_type='log',
                                y_axis_label='fold-change',
                                x_axis_label='IPTG [µM]',
                                title='    O2 | ∆εRA = -13.9 kT    ',
                                y_range=[-0.1, 1.1])
    p_O3 = bokeh.plotting.figure(width=400, height=300, 
                                x_axis_type='log',
                                y_axis_label='fold-change',
                                x_axis_label='IPTG [µM]',
                                title='    O3 | ∆εrA = -9.7 kT    ',
                                y_range=[-0.1, 1.1])

    # ##########################
    # DATA PLOTS 
    # ##########################
    axes = {'O1':p_O1, 'O2':p_O2, 'O3':p_O3}
    rep_colors = {22:pboc['red'], 60:pboc['blue'], 124:pboc['green'],
                  260:pboc['purple'], 1220:pboc['dark_brown'], 1740:pboc['dark_green']}
    for g, d in data.groupby(['operator', 'repressors']):
        # Compute summary statistics
        agg = d.groupby(['IPTGuM']).agg(('mean', 'sem')).reset_index() 

        # Assign axis and plot data 
        _ax = axes[g[0]]
        _ax.circle(agg['IPTGuM'], agg['fold_change']['mean'], color=rep_colors[g[1]])

        # Compute errorbars
        xerrs = [(agg.iloc[i]['IPTGuM'], 
                  agg.iloc[i]['IPTGuM']) for i in range(len(agg))]
        yerrs = [(agg.iloc[i]['fold_change']['mean']-\
                  agg.iloc[i]['fold_change']['sem'], 
                  agg.iloc[i]['fold_change']['mean']+\
                  agg.iloc[i]['fold_change']['sem']) for i in range(len(agg))]
        _ax.multi_line(xerrs, yerrs, line_width=1, color=rep_colors[g[1]])

    lay = layout([[p_O1, p_O2, p_O3]])
    tab = bokeh.models.widgets.Panel(child=lay, title='Razo-Mejia et al. 2018 - Proof of Principle')
    return tab