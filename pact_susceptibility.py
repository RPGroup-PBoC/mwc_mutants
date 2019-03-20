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
bokeh.plotting.output_file("pact_susceptibility.html")

# Define the parameter ranges
n_points = 500

# Define the sliders
ka_slider = Slider(title='Ka [µM]', value=139, start=0.01, end=1000, step=0.01)
ki_slider = Slider(title='Ki [µM]', value=0.53, start=0.01, end=200, step=0.01)
epAI_slider = Slider(title='allosteric energy difference [kT]', value=4.5, 
                    start=-10, end=10, step=0.1)
c_slider = Slider(title='IPTG [µM]', value=50, start=0.01, end=5000, step=0.1)

# Define the range sliders
ka_range = RangeSlider(title='Ka range [µM]', start=0, end=1000, value=(100,500))
ki_range = RangeSlider(title='Ki range [µM]', start=0.001, end=500, value=(0.1,10))
epAI_range = RangeSlider(title='Allosteric energy range [kT]', start=-20, end=20,
                        value=(-5, 15))
c_range = RangeSlider(title='IPTG concentration [µM]', start=0.0001, end=1000,
                       value=(0, 1000))

# Define the source data 
source = ColumnDataSource(data={'ka':[], 'ki':[], 'epAI':[], 'c':[],
                                'dF_depAI':[], 'dF_dka':[], 'dF_dki':[], 
                                'dF_dc':[]})

callback_args = {'source':source, 'nPoints':n_points, 'Ka':ka_slider, 'Ki':ki_slider,
                 'epAI':epAI_slider, 'c':c_slider, 'nSites':2,
                 'rangeKa':ka_range, 'rangeKi':ki_range, 'rangeepAI':epAI_range,
                 'rangeC':c_range}


# Instantiate the figures 
p_epAI = bokeh.plotting.figure(width=425, height=300, 
                               x_axis_label='allosteric energy difference [kT]',
                               y_axis_label='∂∆F/∂∆ε_AI')
p_Ka = bokeh.plotting.figure(width=425, height=300, 
                            x_axis_type='log', x_axis_label = 'Ka [µM]',
                            y_axis_label='∂∆F/∂Ka')
p_Ki = bokeh.plotting.figure(width=425, height=300, 
                            x_axis_type='log', x_axis_label = 'Ki [µM]',
                            y_axis_label='∂∆F/∂Ki')
p_c = bokeh.plotting.figure(width=425, height=300, 
                            x_axis_type='log', x_axis_label = 'IPTG [µM]',
                            y_axis_label='∂∆F/∂c')

# Add initial glyphs. 
p_epAI.line(x='epAI', y='dF_depAI', source=source, color=pboc['red'], line_width=2)
p_Ka.line(x='ka', y='dF_dka', source=source, color=pboc['red'], line_width=2)
p_Ki.line(x='ki', y='dF_dki', source=source, color=pboc['red'], line_width=2)
p_c.line(x='c', y='dF_dc', source=source, color=pboc['red'], line_width=2)

# Define the callback definitions
def callback(cb_list, args=(None,)):
    cb_list.append(" source.change.emit() ")
    return CustomJS(args=args, code=" ".join(cb_list))

variable_definitions=  """
    var data = source.data;
    var nPoints = 1000;

    var ka = Ka.value;
    var ki = Ki.value;
    var epAI = epAI.value;
    var c = c.value;
    var n = nSites.value


    var kaRange = rangeKa.value;
    var kiRange = rangeKi.value;
    var epAIRange = rangeepAI.value;
    var cRange = rangeC.value;

    function dVar(rangeVar) {
        return delta = (rangeVar[1] - rangeVar[0]) / rangeVar.length;
    }
    """

epAI_callback = """
    // Define the functions to compute the derivatives. 
    delta = dVar(epAIRange);
    function derivEpAI(epAIValue) {
        var alloRatio = Math.pow(1 + c / ka, n) / Math.pow(1 + c / ki, n);
        return Math.pow(1 + Math.exp(epAIValue) * alloRatio, -1);
        }

    for (var i = 0; i < nPoints; i++) {
        var deltaVar =  delta * i + epAIRange[0]
        data['epAI'][i] = deltaVar;
        data['dF_depAI'][i] = derivEpAI(deltaVar);
    }
    """
c_callback = """
    var delta = dVar(cRange);
    function derivC(cValue) {
        var numer = 2 * Math.pow(ka, 2) * (ki - ka) * (ki + cValue);
        var denomA = Math.exp(-epAI) * Math.pow(ki, 2) * Math.pow(ka + cValue, 3);
        var denomB = Math.pow(ka, 2) * Math.pow(ki + cValue, 2) * (ka + cValue);
        return numer / (denomA  + denomB);
    }
    for (var i = 0; i < nPoints; i++) {
        var deltaVar = delta * i + cRange[0];
        data['c'][i] = deltaVar;
        data['dF_dc'][i] = derivC(deltaVar);
    }
    """
ka_callback = """
    var  delta = dVar(kaRange);
    function derivKa(kaValue) {
        var numer = -2 * kaValue * c * Math.pow(ki + c, 2);
        var denomA = Math.exp(epAI) * Math.pow(ki, 2) * Math.pow(kaValue + c, 3);
        var denomB = Math.pow(kaValue, 2) * Math.pow(ki + c, 2) * (kaValue + c);
        return numer / (denomA + denomB);
    }

    for (var i = 0; i < nPoints; i++) {
        var deltaVar = delta * i + kaRange[0]
        data['ka'][i] = deltaVar;
        data['dF_dka'][i] = derivKa(deltaVar);
    }
    """
ki_callback = """
    var delta = dVar(kiRange);
    function derivKi(kiValue) {
        var numer = 2 * Math.pow(ka, 2) * c * (kiValue + c);
        var denomA = Math.exp(-epAI) * Math.pow(kiValue, 3) * (ka + c);
        var denomB = kiValue * Math.pow(ka, 2) * Math.pow(kiValue + c, 2);
        return numer / (denomA + denomB);
    }
    for (var i = 0; i < nPoints; i++) {
        var deltaVar = delta * i + kiRange[0];
        data['ki'][i] = deltaVar;
        data['dF_dki'][i] = derivKi(deltaVar); 
    }
    """

master_callback = callback([variable_definitions, epAI_callback,
                            ka_callback, ki_callback, c_callback], args=callback_args)

sliders = [epAI_slider, ka_slider, ki_slider, c_slider]
ranges = [epAI_range, ka_range, ki_range, c_range]

for s, r in zip(sliders, ranges):
    r.callback = master_callback
    s.callback = master_callback

# Generate the layout
constant_parameters = [epAI_slider, ka_slider, ki_slider, c_slider]
range_parameters = [epAI_range, ka_range, ki_range, c_range]
constant_widgets = widgetbox(constant_parameters)
range_widgets = widgetbox(range_parameters)

layout  = bokeh.layouts.layout([[constant_parameters, range_parameters], 
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
