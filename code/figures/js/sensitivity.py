import sys
sys.path.insert(0, '../../../')
import numpy as np
from bokeh.themes import Theme
import mut.viz
import mut.thermo
import bokeh.io
import bokeh.plotting
from bokeh import events
from bokeh.models import ColumnDataSource, Div, CustomJS, CDSView, IndexFilter
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select, Slider, RadioButtonGroup, Button,  RangeSlider
pboc = mut.viz.color_selector('pboc')
bokeh.plotting.output_file("sensitivity.html")


# Define python functions for computing the WT parameter sensitivities
def deriv_epai(ep_value, c, ka, ki, n):
    allo_ratio = ((1 + c/ka) / (1 + c/ki))**n
    return (1 + np.exp(ep_value) * allo_ratio)**-1

def deriv_c(ep, c_value, ka, ki, n):
    numer = n * ((ki * (ka + c_value)) / (ka * (ki + c_value)))**n
    denom_a = (ka + c_value) * (ki + c_value)
    denom_b = ((ki * (ka  c_value) / (ka * (ki + c_value))))**n + np.exp(ep)
    return numer / (denom_a * denom_b)

def deriv_ka(ep, c, ka_value, ki, n):
    numer = c * n * (ki * (ka_value + c) / (ka_value * (ki + c)))**n
    denom_a = ka_value * (ka_value + c)  
    denom_b = ((ki * (ka + c)) / (ka * (ki + c)))**n + np.exp(ep)
    return numer / (denom_a * denom_b)

def deriv_ki(ep, c, ka, ki_value, n):
    numer = -c * n * ((ki_value * (ka + c)) / (ka * (ki_value + c)))**n
    denom_a = ki_value * (ki_value + c) 
    denom_b = ((ki_value * (ka+c)) / (ka * (ki_value + c)))**n + np.exp(ep)
    return numer / (denom_a * denom_b)

# Define the parameter ranges
n_points = 1000

# Don't permit adjustable ranges
ka_range = np.logspace(-2, 3, n_points)
ki_range = np.logspace(-2, 3, n_points)
epAI_range = np.linspace(-10, 10, n_points)
c_range = np.logspace(-3, 3, n_points)

# Define the sliders
ka_slider = Slider(title='Ka [µM]', value=139, start=0.01, end=1000, step=0.01)
ki_slider = Slider(title='Ki [µM]', value=0.53, start=0.01, end=200, step=0.01)
epAI_slider = Slider(title='allosteric energy difference [kT]', value=4.5, 
                    start=-10, end=10, step=0.1)
c_slider = Slider(title='IPTG [µM]', value=50, start=0.01, end=5000, step=0.1)

# Define the range sliders

# Define the source data 
source = ColumnDataSource(data={'ka':ka_range, 'ki':ki_range,
                                'epAI':epAI_range, 'c':c_range,
                                'dF_depAI':np.arange(1, n_points+1, 1), 
                                'dF_dka':np.ones(n_points), 'dF_dki':np.ones(n_points), 
                                'dF_dc':np.ones(n_points)})

# Instantiate the figures 
p_epAI = bokeh.plotting.figure(width=425, height=300, 
                               x_axis_label='allosteric energy difference [kT]',
                               y_axis_label='∂∆F / ∂∆ε_AI',
                               y_range=[0.5, 1.5])
p_Ka = bokeh.plotting.figure(width=425, height=300, 
                            x_axis_type='log', x_axis_label = 'Ka [µM]',
                            y_axis_label='∂∆F / ∂Ka' ,
                            y_range=[-.1, .1])
p_Ki = bokeh.plotting.figure(width=425, height=300, 
                            x_axis_type='log', x_axis_label = 'Ki [µM]',
                            y_axis_label='∂∆F / ∂Ki',
                            y_range=[-1, 1])
p_c = bokeh.plotting.figure(width=425, height=300, 
                            x_axis_type='log', x_axis_label = 'IPTG [µM]',
                            y_axis_label='∂∆F / ∂c',
                            y_range=[-1, 1])

# Add initial glyphs. 
p_epAI.line(x='epAI', y='dF_depAI', source=source, color=pboc['red'], line_width=2)
p_Ka.line(x='ka', y='dF_dka', source=source, color=pboc['red'], line_width=2)
p_Ki.line(x='ki', y='dF_dki', source=source, color=pboc['red'], line_width=2)
p_c.line(x='c', y='dF_dc', source=source, color=pboc['red'], line_width=2)

# Define callback arguments
callback_args = {'source':source, 'nPoints':n_points, 'Ka':ka_slider, 'Ki':ki_slider,
                 'EpAI':epAI_slider, 'C':c_slider, 'nSites':2,
                 'xKa':ka_range, 'xKi':ki_range, 'xepAI':epAI_range,
                 'xC':c_range, 'n_points':n_points, 'nSites':2}


# Define the callback definitions
with open("sensitivity.js") as f:
    callback = CustomJS(args=callback_args, code=f.read())

# Trigger callbacks
sliders = [epAI_slider, ka_slider, ki_slider, c_slider]
plots = [p_epAI, p_Ka, p_Ki, p_c]
for p, s in zip(plots, sliders):
    # p.x_range.range_padding  = 0
    # p.x_range.callback = callback
    s.callback = callback

# Generate the layout
params = [epAI_slider, ka_slider, ki_slider, c_slider]
widgets = widgetbox(params)
# ka_row = bokeh.layouts.column(p_Ka, ka_range)
# epAI_row = bokeh.layouts.column(p_epAI, epAI_range)
# ki_row = bokeh.layouts.column(p_Ki, ki_range)
# c_row = bokeh.layouts.column(p_c, c_range)
layout  = bokeh.layouts.layout([[params], 
                                [p_Ka, p_Ki], [p_epAI, p_c]])

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
