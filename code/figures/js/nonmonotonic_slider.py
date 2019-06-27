# -*- coding: utf-8 -*-
# %%
import numpy as np
import bokeh.io
import bokeh.plotting
import bokeh.layouts
from bokeh.themes import Theme
from bokeh.models import CustomJS, WidgetBox, ColumnDataSource
from bokeh.models.widgets import Slider, RadioButtonGroup
import mut.viz
colors = mut.viz.color_selector('pboc')
bokeh.io.output_notebook()
bokeh.plotting.output_file("nonmonotonicity.html")

# Define the interactions
theta = Slider(title='log\u2081\u2080 θ', start=-3, end=3, step=.01, value=1)
theta1 = Slider(title='log\u2081\u2080 θ\u2081 (modulates Ka)', start=-3, end=3, step=.01, value=1)
theta2 = Slider(title='log\u2081\u2080 θ\u2082 (modulates Ki)', start=-3, end=3, step=.01, value=1)
selector = RadioButtonGroup(labels=['preserve Ka/Ki ratio', 'change Ka/Ki ratio'],
                            active=0)

# Define the data sources
ka = 139
ki = 0.53
ep_ai = 4.5
c_range = np.logspace(-6, 6, 500) 


# Define functions for wild-type and initial cases
def pact(c, ka, ki, ep_ai):
    numer = (1 + c/ka)**2
    denom = numer + np.exp(-ep_ai) * (1 + c/ki)**2
    return numer/denom

def delF(c, ka, ki, ep_ai):
    pact_mut = pact(c, ka, ki, ep_ai)
    pact_wt = pact(c, 139, 0.53, 4.5)
    return - np.log(pact_mut / pact_wt)

def derivF(c, theta1, theta2):
    ka = 139
    ka_mut = theta1 * ka
    ki = 0.53
    ki_mut = theta2 * ki
    ep_ai = 4.5
    wt_numer = ka**2 * (ka - ki) * (c + ki)
    wt_denom = (c + ka) * ((c + ka)**2 * ki**2 +\
                 np.exp(-ep_ai) * ka**2 * (c + ki)**2)
    mut_numer = ka_mut**2 * (ka_mut - ki_mut) * (c + ki_mut)
    mut_denom = (c + ka_mut) * ((c + ka_mut)**2 * ki_mut**2 +\
                 np.exp(-ep_ai) * ka_mut**2 * (c + ki_mut)**2)
    return 2 * np.exp(-ep_ai) * ((mut_numer / mut_denom) - (wt_numer / wt_denom))

wt_pact = pact(c_range, ka, ki, ep_ai)
wt_delF = delF(c_range, ka, ki, ep_ai)
wt_deriv = derivF(c_range, 1, 1)

# Compute the data sources
source = ColumnDataSource({'c_range': c_range, 'x':c_range / ka, 
                            'pact':pact(c_range, 1.1 * ka, 1.1 * ki, ep_ai),
                            'delF':delF(c_range, 1.1 * ka, 1.1 * ki, ep_ai), 
                            'derivF':derivF(c_range, 1.1, 1.1)})


# Define the js callback.
cb = """
// Variable definition
var data = source.data;
var choice = selector.active;
var theta = Math.pow(10, thetaSlider.value);
var theta1 = Math.pow(10, theta1Slider.value);
var theta2 = Math.pow(10, theta2Slider.value);
// Define the wild-type parameters
var ka = 139;
var ki = 0.53;
var ep_ai = 4.5;

// Define the functions
function pact(c, ka, ki, ep_ai) {
    var numer = Math.pow(1 + c / ka, 2);
    var denom = numer + Math.exp(-ep_ai) * Math.pow(1 + c / ki, 2);
    return numer / denom
}

function delF(c, theta1, theta2, ep_ai) {
    var mut_pact = pact(c, theta1 * ka, theta2 * ki, ep_ai);
    var wt_pact = pact(c, ka, ki, ep_ai);
    return Math.log(mut_pact / wt_pact);
}

function derivF(c, theta1, theta2) {
    var ka_mut = theta1 * ka;
    var ki_mut = theta2 * ki;

    wt_numer = Math.pow(ka, 2) * (ka - ki) * (c + ki);
    mut_numer = Math.pow(ka_mut, 2) * (ka_mut - ki_mut) * (c + ki_mut);
    wt_denom = (c + ka) * (Math.pow(c + ka, 2) * 
                Math.pow(ki, 2) + Math.exp(-ep_ai) * 
                Math.pow(ka, 2) * Math.pow(c + ki, 2));
    mut_denom = (c + ka_mut) * (Math.pow(c + ka_mut, 2) * 
                Math.pow(ki_mut, 2) + Math.exp(-ep_ai) * 
                Math.pow(ka_mut, 2) * Math.pow(c + ki_mut, 2))
    return 2 * Math.exp(-ep_ai) * ((mut_numer / mut_denom)  - (wt_numer/wt_denom));
}

// Determine which model to change.
if (choice == 0) {
    var t1 = theta;
    var t2 = theta; 
}

else {
    var t1 = theta1;
    var t2 = theta2;
}

// Update the active probability
for (var i = 0; i < data['c_range'].length; i++) {
data['pact'][i] = pact(data['c_range'][i], t1 * ka, t2 * ki, ep_ai);
data['delF'][i] = delF(data['c_range'][i], t1, t2, ep_ai);
data['derivF'][i] = derivF(data['c_range'][i], t1, t2);}
source.change.emit();
"""

# Define teh callback arguments
cb_args = {'source':source, 'selector':selector, 'thetaSlider':theta,
            'theta1Slider':theta1, 'theta2Slider':theta2}
cb = CustomJS(args=cb_args, code=cb)

# Assign the callback to all widgets
widgets = [selector, theta, theta1, theta2]
for i, w in enumerate(widgets):
    if i == 0:
        w.js_on_change('active', cb)
    else:
        w.js_on_change('value', cb)
    w.callback = cb

# Set up the plot axes
pact_ax = bokeh.plotting.figure(width=300, height=275, x_axis_type='log',
                                x_axis_label='c / Ka', y_axis_label='active probability',
                                y_range=[-0.1, 1.1])

deriv_ax = bokeh.plotting.figure(width=300, height=275, x_axis_type='log',
                                x_axis_label='c / Ka', y_axis_label='\u2202 β∆F/ \u2202 c')

delF_ax = bokeh.plotting.figure(width=600, height=300, x_axis_type='log',
                                x_axis_label='c / Ka', y_axis_label='β∆F',
                                y_range=[-8, 8])

# Plot the wild-type cases 
pact_ax.line(x=c_range / ka, y=wt_pact, color=colors['blue'],
             line_width=2)
delF_ax.line(x=c_range / ka, y=wt_delF, color=colors['blue'], legend='wild type',
            line_width=2)
deriv_ax.line(x=c_range / ka, y=wt_deriv, color=colors['blue'],
            line_width=2)

# Plot the initial mutant cases
pact_ax.line(x='x', y='pact', color=colors['red'], line_width=2, source=source)
deriv_ax.line(x='x', y='derivF', color=colors['red'], line_width=2, source=source)
delF_ax.line(x='x', y='delF', color=colors['red'], line_width=2, source=source,
            legend='mutant')

col1 = bokeh.layouts.column(theta1, theta2)
row1 = bokeh.layouts.row(theta, col1)
ax_row1 = bokeh.layouts.row(pact_ax, deriv_ax)
interactions = bokeh.layouts.column(selector, row1)
plots = bokeh.layouts.column(ax_row1, delF_ax)
fig = bokeh.layouts.column(interactions, plots)

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
bokeh.io.save(fig)

#%%
