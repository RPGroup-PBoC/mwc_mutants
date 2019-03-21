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

# Define the parameter ranges
n_points = 500

# Define the sliders
ka_slider = Slider(title='Ka [µM]', value=139, start=0.01, end=1000, step=0.01)
ki_slider = Slider(title='Ki [µM]', value=0.53, start=0.01, end=200, step=0.01)
epAI_slider = Slider(title='allosteric energy difference [kT]', value=4.5, 
                    start=-10, end=10, step=0.1)
c_slider = Slider(title='IPTG [µM]', value=50, start=0.01, end=5000, step=0.1)

# Define the source data 
source = ColumnDataSource(data={'ka':np.ones(n_points), 'ki':np.ones(n_points), 
                                'epAI':np.arange(0, n_points, 1), 'c':np.ones(n_points),
                                'dF_depAI':np.arange(0, n_points, 1), 
                                'dF_dka':np.ones(n_points), 'dF_dki':np.ones(n_points), 
                                'dF_dc':np.ones(n_points)})

# Instantiate the figures 
p_epAI = bokeh.plotting.figure(width=425, height=300, 
                               x_axis_label='allosteric energy difference [kT]',
                               y_axis_label='∂∆F / ∂∆ε_AI')
p_Ka = bokeh.plotting.figure(width=425, height=300, 
                            x_axis_type='log', x_axis_label = 'Ka [µM]',
                            y_axis_label='∂∆F / ∂Ka')
p_Ki = bokeh.plotting.figure(width=425, height=300, 
                            x_axis_type='log', x_axis_label = 'Ki [µM]',
                            y_axis_label='∂∆F / ∂Ki')
p_c = bokeh.plotting.figure(width=425, height=300, 
                            x_axis_type='log', x_axis_label = 'IPTG [µM]',
                            y_axis_label='∂∆F / ∂c')

# Add initial glyphs. 
p_epAI.line(x='epAI', y='dF_depAI', source=source, color=pboc['red'], line_width=2)
p_Ka.line(x='ka', y='ka', source=source, color=pboc['red'], line_width=2)
p_Ki.line(x='ki', y='ki', source=source, color=pboc['red'], line_width=2)
p_c.line(x='c', y='c', source=source, color=pboc['red'], line_width=2)

# Define callback arguments
callback_args = {'source':source, 'nPoints':n_points, 'ka':ka_slider, 'ki':ki_slider,
                 'epAI':epAI_slider, 'c':c_slider, 'nSites':2,
                 'xKa':p_Ka.x_range, 'xKi':p_Ki.x_range, 'xepAI':p_epAI.x_range,
                 'xC':p_c.x_range, 'n_points':n_points}


# Define the callback definitions
def callback_gen(cb_list, args=callback_args):
   return CustomJS(args=args, code=" ".join(cb_list)) 


fn_callback = """
    // Define the parameter values
    var data = source.data;
    var ka = ka.value;
    var ki = ki.value;
    var epAI = epAI.value;
    var c = c.value;
    var n = nSites.value;

    function linspace(start, stop, n) {
        var x = [];
        var currValue = start;
        var step = (stop - start) / (n - 1);
        for (var i = 0; i < n; i++) {
            x.push(currValue);
            currValue += step;
        }
        return x;
    }

    function logspace(start, stop, n) {
        var x = [];
        var currValue = start;
        var step = (stop - start) / n;
        for (var i = 0; i < n; i++) {
            x.push(Math.pow(10, currValue));
            currValue += step   
        }
    }

    function derivEpAI(epValue, c, ka, ki, n){
        var allo_ratio =  Math.pow(1 + c / ka, n) / Math.pow(1 + c/ki, n);
        return Math.pow(1 + Math.exp(epValue) * allo_ratio, -1);
    }

    function derivC(ep, cValue,  ka, ki, n){
        var numer = n * Math.pow((ki * (ka + cValue))/(ka * (ki + cValue)),n) * (ka - ki);
        var denom_a = (ka + cValue) * (ki + cValue);
        var denom_b = (Math.pow((ki * (ka + cValue) / (ka * (ki + cValue))), n) + Math.exp(ep));
        return numer / (denom_a * denom_b);
    }

    function derivKa(ep, c, kaValue, ki, n) {
        var numer = c * n * Math.pow((ki * (kaValue + c))/ (kaValue * (ki + c)), n);
        var denom_a = kaValue * (kaValue + c)
        var denom_b = Math.pow((ki * (ka + c))/ (ka * (ki + c)), n) + Math.exp(ep);
        return numer / (denom_a * denom_b);
    }

    function derivKi(ep, c, ka, kiValue, n) { 
        var numer = -c * n * Math.pow((kiValue * (ka + c))/ (ka * (kiValue + c)), n);
        var denom_a = kiValue * (kiValue + c)
        var denom_b = Math.pow((kiValue * (ka + c)) / (ka * (kiValue + c)), n)  + Math.exp(ep)
        return numer / (denom_a * denom_b);
    }

    """
dF_depAI_callback = """
    data['epAI'] = linspace(xepAI.start, xepAI.end, n_points); 
    for (var i = 0; i < n_points; i++) {
   data['dF_depAI'][i] = derivEpAI(data['epAI'][i], c, ka, ki, n);}
   source.change.emit()
    """
df_dKa_callback = """
    data['ka'] = linspace(xKa.start, xKa.end, n_points);
    """
# with open("./sensitivity.js") as f:
    # callback_str = f.read()
# callback = CustomJS(args=callback_args, code=callback_str)

# Trigger callbacks
sliders = [epAI_slider, ka_slider, ki_slider, c_slider]
plots = [p_epAI, p_Ka, p_Ki, p_epAI]
for p, s in zip(plots, sliders):
    p.x_range.callback = callback_gen([fn_callback, dF_depAI_callback])
    s.js_on_change('value', callback_gen([fn_callback, dF_depAI_callback]))

# Generate the layout
params = [epAI_slider, ka_slider, ki_slider, c_slider]
widgets = widgetbox(params)

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
