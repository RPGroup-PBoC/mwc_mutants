# -*- coding: utf-8 -*-
# %%
import sys
sys.path.insert(0, '../../../')
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
from bokeh.models.widgets import TextInput, Slider, Toggle, Button, RadioButtonGroup
from bokeh.models import CustomJS, ColumnDataSource
import bokeh.transform
import bokeh.palettes
import bokeh.layouts
import mut.thermo
constants = mut.thermo.load_constants()
bokeh.plotting.output_file("param_selector.html")

# Load the data set (Razo-Mejia et al 2018)
data = pd.read_csv('../../../data/csv/RazoMejia2018_data.csv')

# Clean the data
data = data[data['repressors'] > 0]
data['repressors'] *= 2
data['repressors'] = data['repressors'].astype('str')
data.rename(columns={'fold_change_A':'fold_change', 'IPTG_uM':'IPTGuM'}, inplace=True)

# Assign the color idx
pal = bokeh.transform.factor_cmap('repressors', palette=bokeh.palettes.Dark2_6, 
                    factors=data['repressors'].unique())

# Define the raw data sets as sources
O1_data = data[data['operator']=='O1']
O2_data = data[data['operator']=='O2']
O3_data = data[data['operator']=='O3']

# Define the raw data sets as sources
O1_summ = O1_data.groupby(['repressors', 
                    'IPTGuM'])['fold_change'].agg(('mean', 'sem')).reset_index()
O2_summ = O2_data.groupby(['repressors', 
                    'IPTGuM'])['fold_change'].agg(('mean', 'sem')).reset_index()
O3_summ = O3_data.groupby(['repressors', 
                    'IPTGuM'])['fold_change'].agg(('mean', 'sem')).reset_index()
# Compute the errorbars
O1_summ['low'] = O1_summ['mean'] - O1_summ['sem']
O1_summ['high'] = O1_summ['mean'] + O1_summ['sem']
O2_summ['low'] = O2_summ['mean'] - O2_summ['sem']
O2_summ['high'] = O2_summ['mean'] + O2_summ['sem']
O3_summ['low'] = O3_summ['mean'] - O3_summ['sem']
O3_summ['high'] = O3_summ['mean'] + O3_summ['sem']

# Generate the sources
O1_summ = ColumnDataSource(O1_summ)
O2_summ = ColumnDataSource(O2_summ)
O3_summ = ColumnDataSource(O3_summ)

# Generate theory sources with initial states 
c_range = np.logspace(-2, 4, 100)
ka_init = 139
ki_init = 0.53
krr_init = np.exp(-4.5)
sources = []
for o in data['operator'].unique():
    xs = []
    ys = []
    cs = []
    for r in data['repressors'].unique():
        fc = mut.thermo.SimpleRepression(float(r), ep_r=constants[o], ka=ka_init, ki=ki_init,
                                        ep_ai=4.5, effector_conc=c_range).fold_change()
        xs.append(c_range)
        ys.append(fc)
        cs.append(r)
    sources.append({'IPTGuM': xs, 'fold_change':ys, 'repressors':cs})

O1_theo = ColumnDataSource(sources[1])
O2_theo = ColumnDataSource(sources[0])
O3_theo = ColumnDataSource(sources[2])

#%%
# Instantiate the direct text inputs
ka_input = TextInput(value="139", title='Ka [µM]')
ki_input = TextInput(value="0.53", title='Ki [µM]')
krr_input = TextInput(value="0.001", title='K_RR* (exp{-βΔε})')

# Instantiate the slider inputs
ka_slider = Slider(start=0, end=500, step=0.01, value=139, title='Ka [µM]')
ki_slider = Slider(start=0, end=500, step=0.01, value=0.53, title='Ki [µM]')
krr_slider = Slider(start=-10, end=10, step=0.01, value=0.01, title='K_RR* (exp{-βΔε})')

# Instantiate the selector buttons
selector = RadioButtonGroup(labels=['Use Numeric Input', 'Use Slider Input'], active=0)

# Execute button
execute = Button(label='Update Plots', button_type='success')

# Generate layouts for widgets
box1 = bokeh.layouts.widgetbox([ka_input, ki_input,
                                 krr_input])
box2 = bokeh.layouts.widgetbox([ka_slider, ki_slider, 
                                krr_slider])
row1 = bokeh.layouts.row(box1, box2)
widget_layout = bokeh.layouts.column(selector, row1)


# load the callbacks
with open('param_selector.js') as file:
    cb = file.read()

callback = CustomJS(args=dict(o1_source=O1_theo, o2_source=O2_theo, 
                            o3_source=O3_theo, Ka_numer=ka_input, 
                            Ki_numer=ki_input, Krr_numer=krr_input,
                            Ka_slider=ka_slider, Ki_slider=ki_slider, 
                            Krr_slider=krr_slider, Radio=selector, c_range=c_range), code=cb)

# Apply the callback to the update plots
ka_input.callback = callback
ki_input.callback = callback
krr_input.callback = callback
ka_input.js_on_change('value', callback)
ki_input.js_on_change('value', callback)
krr_input.js_on_change('value', callback)
ka_slider.callback = callback
ki_slider.callback = callback
krr_slider.callback = callback
ka_slider.js_on_change('value', callback)
ki_slider.js_on_change('value', callback)
krr_slider.js_on_change('value', callback)

# Set up the figure canvas
o1_ax = bokeh.plotting.figure(width=300, height=300, x_axis_type='log',
                            x_axis_label='IPTG [µM]', y_axis_label='fold-change',
                            title='Operator O1 : Δε_RA = -15.3 kT')
o2_ax = bokeh.plotting.figure(width=300, height=300, x_axis_type='log',
                            x_axis_label='IPTG [µM]', y_axis_label='fold-change',
                            title='Operator O2 : Δε_RA = -13.9 kT')
o3_ax = bokeh.plotting.figure(width=300, height=300, x_axis_type='log',
                            x_axis_label='IPTG [µM]', y_axis_label='fold-change',
                            title='Operator O3 : Δε_RA = -9.7 kT')

# Plot the theory
o1_ax.multi_line(xs='IPTGuM', ys='fold_change', line_width=1, source=O1_theo, line_color=pal)
o2_ax.multi_line(xs='IPTGuM', ys='fold_change', line_width=1, source=O2_theo, line_color=pal)
o3_ax.multi_line(xs='IPTGuM', ys='fold_change', line_width=1, source=O3_theo, line_color=pal)

# Plot the summarized data
o1_ax.segment(x0='IPTGuM', x1='IPTGuM', y0='low', y1='high', line_width=1.5, 
                source=O1_summ, legend='repressors', line_color=pal)
o1_ax.circle(x='IPTGuM', y='mean', source=O1_summ, legend='repressors', line_color=pal, fill_color='white', size=5)
o2_ax.segment(x0='IPTGuM', x1='IPTGuM', y0='low', y1='high', line_width=1.5, 
                source=O2_summ, line_color=pal)
o2_ax.circle(x='IPTGuM', y='mean', source=O2_summ, fill_color='white', line_color=pal, size=5)
o3_ax.segment(x0='IPTGuM', x1='IPTGuM', y0='low', y1='high', line_width=1.5, 
                source=O3_summ, line_color=pal)
o3_ax.circle(x='IPTGuM', y='mean', source=O3_summ, fill_color='white', line_color=pal, size=5)

# Format the legends
o1_ax.legend.location = 'top_left'
o1_ax.legend.title = 'R'

# Define the plot layouts
plots = bokeh.layouts.row(o1_ax, o2_ax, o3_ax)
layout = bokeh.layouts.column(widget_layout, plots)
bokeh.io.save(layout)
#%%


#%%
