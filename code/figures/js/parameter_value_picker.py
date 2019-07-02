# -*- coding: utf-8 -*-
"""
This scripts generates an interactive figure to explore how different parameter
values control the shapes of wild-type LacI induction curves.
"""
# %%
# ##############################################################################
# IMPORTS AND CONSTANT DEFINITIONS
# ##############################################################################
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
from bokeh.models.widgets import (TextInput, Slider, Toggle, 
                                 Button, RadioButtonGroup, Dropdown, Div)
from bokeh.models import CustomJS, ColumnDataSource
import bokeh.transform
import bokeh.palettes
import bokeh.layouts
from bokeh.themes import Theme
import mut.thermo
constants = mut.thermo.load_constants()
constants['Oid'] = -17
bokeh.plotting.output_file("param_selector.html")
# bokeh.io.output_notebook()

# ##############################################################################
# DATA LOADING, CLEANING, AND SOURCE DEFINITION
# ##############################################################################
# Load the data set (Razo-Mejia et al 2018) and clean
data = pd.read_csv('../../../data/csv/RazoMejia2018_data.csv')
data = data.sort_values(['repressors'])
data = data[data['repressors'] > 0]
data['repressors'] *= 2
data['repressors'] = data['repressors'].astype('str')
data.rename(columns={'fold_change_A':'fold_change', 'IPTG_uM':'IPTGuM'}, inplace=True)

# Load the fugacity data and leakiness data
old_gods = pd.read_csv('../../../data/csv/Garcia2011_Brewster2014_data.csv')
garcia = old_gods[old_gods['author']=='garcia']
brewster = old_gods[old_gods['author']=='brewster']
fug_data = pd.read_csv('../../../data/csv/fugacity_data.csv')
fug_data = fug_data[fug_data['operator']=='O1']
fug_data['N'] = fug_data['N'].astype(str)

# Load the global sampling statistics for three different choices of KRR
razo_stats = pd.read_csv('../../../data/csv/razomejia_epai_summary.csv')

# Assign the color idx
pal = bokeh.transform.factor_cmap('repressors', palette=bokeh.palettes.Dark2_6, 
                    factors=data['repressors'].unique())

# Define the raw data sets as sources
O1_data = data[data['operator']=='O1']
O2_data = data[data['operator']=='O2']
O3_data = data[data['operator']=='O3']
garcia_source = ColumnDataSource(garcia)
brewster_source = ColumnDataSource(garcia)
fug_source = ColumnDataSource(fug_data)

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
for o in ['O1', 'O2', 'O3', 'Oid']:
    xs = []
    ys = []
    cs = []
    for r in data['repressors'].unique():
        fc = mut.thermo.SimpleRepression(float(r), ep_r=constants[o], ka=ka_init, ki=ki_init,
                                        ep_ai=4.5, effector_conc=c_range).fold_change()
        xs.append(c_range)
        ys.append(fc)
        cs.append([r])
    sources.append({'IPTGuM': xs, 'fold_change':ys, 'repressors':cs})
O1_theo = ColumnDataSource(sources[1])
O2_theo = ColumnDataSource(sources[0])
O3_theo = ColumnDataSource(sources[2])

# Leakiness plots
rep_range = np.logspace(0, 4, 100)
xs = []
ys = []
cs = []
for o in ['O1', 'O2', 'O3', 'Oid']:
    xs.append(rep_range)
    pact = (1 + np.exp(-constants['ep_AI']))
    Ra = pact * rep_range / 4.6E6
    ys.append((1 + Ra * np.exp(-constants[o]))**-1)
    cs.append(o)
leak_theo = ColumnDataSource({'repressors':xs, 'leakiness':ys, 'operator':cs})

pal2 = bokeh.transform.factor_cmap('repressors', palette=bokeh.palettes.Dark2_6, 
                                   factors=data.repressors.unique())

# Compute the fugacity data
xs = []
ys = []
cs = []
for n in [64, 52]:
    xs.append(rep_range)
    x = np.exp(15.3)
    pact = 1 / (1 + np.exp(-constants['ep_AI']))
    r = pact * rep_range
    B = rep_range * (1 + x) - n * x - 4.6E6
    A = x * (r - n - 4.6E6)
    numer = -B - np.sqrt(B * B - 4 * A * r)
    denom = 2 * A
    ys.append(1 / (1 + (numer / denom) * np.exp(15.3)))
    cs.append(str(n))
fug_theo = ColumnDataSource({'repressors':xs, 'leakiness':ys, 'N':cs})

#%%
# Instantiate the direct text inputs
ka_input = TextInput(value="139", title='Ka [µM]')
ki_input = TextInput(value="0.53", title='Ki [µM]')
krr_input = TextInput(value="0.001", title='K_RR* (exp{-βΔε})')

# Instantiate the slider inputs
ka_slider = Slider(start=0, end=500, step=0.01, value=139, title='Ka [µM]')
ki_slider = Slider(start=0, end=500, step=0.01, value=0.53, title='Ki [µM]')
krr_slider = Slider(start=0.000001, end=10, step=0.01, value=0.01, title='K_RR* (exp{-βΔε})')

# Instantiate the selector buttons
selector = RadioButtonGroup(labels=['Use Numeric Input', 'Use Slider Input', 'Use Literature Values'], active=0)
drop = Dropdown(label='Select Literature Source', menu=[('Razo-Mejia et al. 2018', 'razo'), 
                                                        ('Daber, Sharp, and Lewis 2009 (Δε_RI = 0)', 'one_not_enough_x0'),
                                                        ('Daber, Sharp, and Lewis 2009 (Δε_RI = Δε_RA / 100,000)', 'one_not_enough_x100000'),
                                                        ('Daber, Sharp, and Lewis 2009 (Δε_RI = Δε_RA / 10,000)', 'one_not_enough_x10000'),
                                                        ('Daber, Sharp, and Lewis 2009 (Δε_RI = Δε_RA / 1000)', 'one_not_enough_x1000'),
                                                        ('Daber, Sharp, and Lewis 2009 (Δε_RI = Δε_RA / 100)', 'one_not_enough_x100'),
                                                        ('Daber, Sochor, and Lewis 2011', 'daber_muts'),
                                                        ("O'Gorman et al. 1980", 'ogorman')],
                                                        button_type='warning')

# Describe the parameter values.
vals = Div(text="""Literature values appear here when selection is made.""", width =300, height=150)
div1 = Div(text="""
<b>This interactive figure allows you to explore how changing parameter values of the simple repression thermodynamic model match <i>in vivo</i> data collected in the Phillips Group. To use, first select the mode of entry (Numeric, Slider, or from Literature Values) followed by your selection of parameters. The data should update as soon as a selection is made.</b>
""", width = 900, height=50)
div2 = Div(text=""" <center><b>These data come from Razo-Mejia <i>et al.</i> 2018 in which
the fold-change in gene expression was measured <i> in vivo</i>.</b></center>""", width=900, height=75)
div3 = Div(text=""" <center><b>These data come Garcia and Phillips 2011 and Brewster <i>et al.</i> 2014. They show the fold-change in gene expression measured for LacI tetramers (circles) and dimers (diamonds) for four different operators in the absence of inducer. In these experiments, there is a single specific binding site in the genome. the fold-change in gene expression was measured <i> in vivo</i>.</b></center>""", width=450, height=75)
div4 = Div(text="""<center><b>These data come from Brewster <i>et al.</i> 2014 in which the fold-chang in gene expression for a simple repression architecture with the native O1 operator sequence was present in the cell with either 64 (green) or 52 (orange) copies of the promoter.</b></center>""", width=400, height=75)

# Execute button
execute = Button(label='Update Plots', button_type='success')

# Generate layouts for widgets
box1 = bokeh.layouts.widgetbox([ka_input, ki_input,
                                 krr_input])
box2 = bokeh.layouts.widgetbox([ka_slider, ki_slider, 
                                krr_slider])

box3 = bokeh.layouts.widgetbox([drop, vals])
row1 = bokeh.layouts.row(box1, box2, box3)
widget_layout = bokeh.layouts.column(selector, row1)


# load the callbacks
with open('param_selector.js') as file:
    cb = file.read()

callback = CustomJS(args=dict(o1_source=O1_theo, o2_source=O2_theo, 
                            o3_source=O3_theo, Ka_numer=ka_input, 
                            Ki_numer=ki_input, Krr_numer=krr_input,
                            Ka_slider=ka_slider, Ki_slider=ki_slider, 
                            Krr_slider=krr_slider, Radio=selector, 
                            c_range=c_range, Drop=drop, Desc=vals, 
                            rep_range=rep_range, leak_source=leak_theo,
                            fug_source=fug_theo), code=cb)

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
drop.js_on_change("value", callback)

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
old_ax = bokeh.plotting.figure(width=450, height=300, 
                                  x_axis_type='log', y_axis_type='log',
                                  x_axis_label='repressors per cell', y_axis_label='leakiness',
                                  title="Single promoter, no inducer")
fug_ax = bokeh.plotting.figure(width=450, height=300, 
                                  x_axis_type='log', y_axis_type='log',
                                  x_axis_label='repressors per cell', y_axis_label='leakiness',
                                  title='variable promoter copy number (O1 operator)',
                                  x_range = [1, 1E3])

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

# Plot the brewster and garcia data
pal3 = bokeh.transform.factor_cmap(field_name='operator', palette=bokeh.palettes.Dark2_6, factors=['O1', 'O2', 'O3', 'Oid'])
pal4 = bokeh.transform.factor_cmap(field_name='N', palette=bokeh.palettes.Dark2_6, factors=fug_data['N'].unique())
old_ax.circle(x='repressor', y='fold_change', fill_color='white', 
            legend='operator', line_color=pal3, size=6, source=garcia_source)
old_ax.diamond(x='repressor', y='fold_change', fill_color='white', line_color=pal3, size=8, source=brewster)
fug_ax.circle(x='repressor', y='fold_change', source=fug_source, fill_color='white', 
             line_color=pal4, size=6, legend='N')

# Plot the Leakiness and fugacity theory
old_ax.multi_line(xs='repressors', ys='leakiness', color=pal3, source=leak_theo, line_width=1)
fug_ax.multi_line(xs='repressors', ys='leakiness', color=pal4, source=fug_theo, line_width=1)

# Format the legends
o1_ax.legend.location = 'top_left'
o1_ax.legend.title = 'R'
old_ax.legend.location = 'bottom_left'
fug_ax.legend.title = 'number of promoters'
fug_ax.legend.location = 'bottom_left'

# Define the plot layouts
razo_plots = bokeh.layouts.row(o1_ax, o2_ax, o3_ax)
old_god_plots = bokeh.layouts.row(old_ax, fug_ax)
div_row = bokeh.layouts.row(div3, div4)
layout = bokeh.layouts.column(div1, widget_layout, razo_plots, old_god_plots)
# bokeh.io.show(layout)

# #############################
# THEME DETAILS
# ############################
theme_json = {'attrs':
            {'Figure': {
                'background_fill_color': '#E3DCD0',
                'outline_line_color': '#FFFFFF',
            },
            'Axis': {
            'axis_line_color': "white",
            'major_tick_in': 7,
            'major_tick_line_width': 2.5,
            'major_tick_line_color': "white",
            'minor_tick_line_color': "white",
            'axis_label_text_font': 'Helvetica',
            'axis_label_text_font_style': 'normal'
            },
            'Grid': {
                'grid_line_color': None,
            },
            'Legend': {
                'background_fill_color': '#E3DCD0',
                'border_line_color': '#FFFFFF',
                'border_line_width': 1.5,
                'background_fill_alpha': 0.5
            },
            'Text': {
                'text_font_style': 'normal',
               'text_font': 'Helvetica'
            },
            'Title': {
                'background_fill_color': '#FFEDC0',
                'text_font_style': 'normal',
                'align': 'center',
                'text_font': 'Helvetica',
                'offset': 2,
            }}}

theme = Theme(json=theme_json)
bokeh.io.curdoc().theme = theme
bokeh.io.save(layout)
# bokeh.io.show(layout)


#%%
