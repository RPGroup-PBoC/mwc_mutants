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
    allo_ratio = ((ka * (ki + c))/ (ki * (ka + c)))**n
    return allo_ratio / (allo_ratio + np.exp(ep_value))

def deriv_c(ep, c_value, ka, ki, n):
    numer = -n * ((ka * (ki + c_value)) / (ki * (ka + c_value)))**n * (ka - ki)
    denom_a = (ka + c_value) * (ki + c_value)
    denom_b = ((ka * (ki + c_value) / (ki * (ka + c_value))))**n + np.exp(ep)
    return numer / (denom_a * denom_b)

def deriv_ka(ep, c, ka_value, ki, n):
    numer = -c * n * (ka_value * (ki + c) / (ki * (ka_value + c)))**n
    denom_a = ka_value * (ka_value + c)  
    denom_b = ((ka_value * (ki + c)) / (ki * (ka_value + c)))**n + np.exp(ep)
    return numer / (denom_a * denom_b)

def deriv_ki(ep, c, ka, ki_value, n):
    numer = c * n * ((ka * (ki_value + c)) / (ki_value * (ka + c)))**n
    denom_a = ki_value * (ki_value + c) 
    denom_b = ((ka * (ki_value + c)) / (ki_value * (ka + c)))**n + np.exp(ep)
    return numer / (denom_a * denom_b)

# Don't permit adjustable ranges
n_points = 300
log_ka_range = np.linspace(-4, 4, n_points)
log_ki_range = np.linspace(-4, 4, n_points)
epAI_range = np.linspace(-8, 8, n_points)
log_c_range = np.linspace(-4, 4, n_points)

# Define the sliders
ka_slider = Slider(title='log10(Ka / 1 µM)', value=np.log10(139), start=log_ka_range[0], 
                                     end=log_ka_range[-1], step=0.01, 
                                     bar_color=pboc['light_red'])
ki_slider = Slider(title='log10(Ki / 1µM)', value=np.log10(0.53), start=log_ki_range[0], 
                   end=log_ki_range[-1], step=0.01, bar_color=pboc['light_red'])
epAI_slider = Slider(title='allosteric energy difference [kT]', value=4.5, 
                    start=-10, end=10, step=0.1, bar_color=pboc['light_red'])
c_slider = Slider(title='log10(IPTG / 1µM)', value=np.log10(50), start=log_c_range[0], 
                    end=log_c_range[-1], step=0.01, bar_color=pboc['light_red'])

# # Define the source data 
source = ColumnDataSource(data={'ka':10**log_ka_range, 'ki':10**log_ki_range,
                                'epAI':epAI_range, 'c':10**log_c_range,
                                'dF_depAI':deriv_epai(epAI_range, 50, 139, 0.53, 2), 
                                'dF_dka': deriv_ka(4.5, 50, 10**log_ka_range, 0.53, 2),
                                'dF_dki':deriv_ki(4.5, 50, 139, 10**log_ki_range, 2), 
                                'dF_dc':deriv_c(4.5, 10**log_c_range, 139, 0.53, 2)})

# Instantiate the figures 
p_epAI = bokeh.plotting.figure(width=425, height=300, 
                               x_axis_label='allosteric energy difference [kT]',
                               y_axis_label='∂∆F / ∂∆ε_AI',
                               y_range=[-0.1, 1.1])
p_Ka = bokeh.plotting.figure(width=425, height=300, 
                            x_axis_type='log',
                            x_axis_label = 'Ka [µM]',
                            y_axis_label='∂∆F / ∂Ka' )
p_Ki = bokeh.plotting.figure(width=425, height=300, 
                            x_axis_type='log',
                            x_axis_label = 'Ki [µM]',
                            y_axis_label='∂∆F / ∂Ki')
p_c = bokeh.plotting.figure(width=425, height=300, 
                            x_axis_type='log',
                            x_axis_label = 'IPTG [µM]',
                            y_axis_label='∂∆F / ∂c')
                            

# Add initial glyphs. 
p_epAI.line(x='epAI', y='dF_depAI', source=source, color=pboc['red'], line_width=2)
p_Ka.line(x='ka', y='dF_dka', source=source, color=pboc['red'], line_width=2)
p_Ki.line(x='ki', y='dF_dki', source=source, color=pboc['red'], line_width=2)
p_c.line(x='c', y='dF_dc', source=source, color=pboc['red'], line_width=2)

# Define callback arguments
callback_args = {'source':source, 'Ka':ka_slider, 'Ki':ki_slider,
                 'EpAI':epAI_slider, 'C':c_slider, 'nSites':2,
                 'nPoints':n_points}

# Load the callback files
with open("sensitivity.js") as f:
    sensitivity = f.read()

general_callback = CustomJS(args=callback_args, code=sensitivity)
_reset_callback = CustomJS(args=callback_args, code="""
    var ka = 139; 
    var ki = 0.53;
    var c = 50;
    var epAI = 4.5;
    var n = 2;
    Ka.value = Math.log10(ka);
    Ki.value = Math.log10(ki);
    EpAI.value = epAI;
    C.value = Math.log10(c);
 """)

# Define the reset button
reset = Button(label='Reset to default', callback=_reset_callback)
reset.js_on_click(general_callback)

# Trigger callbacks
sliders = [epAI_slider, ka_slider, ki_slider, c_slider]
for s in sliders:
    s.callback =  general_callback

# Generate the layout
params = [epAI_slider, ka_slider, ki_slider, c_slider]
layout  = bokeh.layouts.layout([[reset], [params[:2], params[2:]], 
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