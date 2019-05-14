import numpy as np
import pandas as pd
import mut.thermo
import mut.viz
import bokeh.io
import bokeh.plotting
import bokeh.models
from bokeh import events
from bokeh.models import ColumnDataSource, Div, CustomJS, CDSView, IndexFilter
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select, Slider, RadioButtonGroup, Button
constants = mut.thermo.load_constants()
colors = mut.viz.color_selector('pboc')

# Load mitch's data
daber_data = pd.read_csv('../../../data/csv/Daber2011_data.csv')
wt = daber_data[daber_data['mutant']=='wt']
ind_data = daber_data[(daber_data['mutant']=='F161T') | (daber_data['mutant']=='Q291V') |
                 (daber_data['mutant']=='Q291R') | (daber_data['mutant']=='Q291K')] 

# Define the javascript
js = """
// Define the constants
var data = source.data;
var c = data['c'];
var R = R_slider.value
var M = M_slider.value
var ka = 139;
var ki = 0.53
var ep_ai = 4.5
var n = 2;
var Nns = 46000000;

// Hard-code the mutant Ks
var F161T_ka = 300;
var F161T_ki = 13;
var F161T_epAI = -1;
var Q291V_ka = 1100;
var Q291V_ki = 53;
var Q291V_epAI = 0;
var Q291R_ka = 135;
var Q291R_ki = 143;
var Q291R_epAI = -2.4;
var Q291K_ka = 1450;
var Q291K_ki = 330;
var Q291K_epAI = -3.16;

// Define function to compute the fugacity. 
function computePact(c, ka, ki, ep_ai, n) {
    var numer = Math.pow(1 + c / ki, n);
    var denom = Math.pow(1 + c / ka, n);
    var prob = Math.pow(1 + Math.exp(-ep_ai) * numer / denom, -1);
    return prob;
    }
function computeFoldChange(R, M, pact) {
    var r = R * pact;
    var x = Math.exp(15.3); // Defined value for O1
    var B = (r * (1 + x) - M * x - Nns) 
    var C_term = r;
    var A = x * (r - M - Nns);
    var numer = -B - Math.sqrt(B*B - 4 * A * C_term);
    var denom = 2 * A;
    var fugacity = numer / denom;
    return Math.pow(1 + fugacity * x, -1); 
    }
    
// Evaluate the fold-change. 
var fc = [];
var F161T = [];
var Q291V = [];
var Q291R = [];
var Q291K = [];
for (var i = 0; i < c.length; i++) {
   var pact = computePact(c[i], ka, ki, ep_ai, n); 
   fc[i] = computeFoldChange(R, M, pact);
 
   F161T_pact = computePact(c[i], F161T_ka, F161T_ki, F161T_epAI, n);
   Q291V_pact = computePact(c[i], Q291V_ka, Q291V_ki, Q291V_epAI, n);
   Q291R_pact = computePact(c[i], Q291R_ka, Q291R_ki, Q291R_epAI, n);
   Q291K_pact = computePact(c[i], Q291K_ka, Q291K_ki, Q291K_epAI, n); 
   F161T[i] = computeFoldChange(R, M, F161T_pact);
   Q291V[i] = computeFoldChange(R, M, Q291V_pact); 
   Q291R[i] = computeFoldChange(R, M, Q291R_pact); 
   Q291K[i] = computeFoldChange(R, M, Q291K_pact);    
   } 


source.data['fold_change'] = fc;
source.data['F161T'] = F161T;
source.data['Q291V'] = Q291V;
source.data['Q291R'] = Q291R;
source.data['Q291K'] = Q291K;
source.change.emit()
"""

# Define the plotting backbone
bokeh.plotting.output_file("fugacity_explorer.html")
c_range = np.logspace(-2,6, 500)
p = bokeh.plotting.figure(width=600, height=400, x_axis_type='log',
                         x_axis_label='IPTG [µM]', y_axis_label='fold-change')
p2 = bokeh.plotting.figure(width=600, height=400, x_axis_type='log',
                         x_axis_label='IPTG [µM]', y_axis_label='fold-change')

# Set up the column data source for the slider.
ones = np.ones(len(c_range));
source = ColumnDataSource(data={'c':c_range, 'fold_change':ones,
                                'F161T':ones, 'Q291V':ones, 'Q291R':ones, 
                                'Q291K':ones})
 
# Add various glyphs
p.circle(wt['IPTGuM'], wt['fold_change'], legend='Wild Type Data', size=8)
p.line(x='c', y='fold_change', source=source, line_width=4, legend='Theory')

colors = {'F161T':colors['red'], 'Q291V':colors['green'], 'Q291K':colors['purple'], 'Q291R':colors['dark_brown']}
for g, d in ind_data.groupby('mutant'):
    p2.circle(d['IPTGuM'], d['fold_change'], color=colors[g], size=8, legend=g)

p2.line(x='c', y='F161T', source=source, color=colors['F161T'], line_width=4)
p2.line(x='c', y='Q291V', source=source, color=colors['Q291V'], line_width=4)
p2.line(x='c', y='Q291R', source=source, color=colors['Q291R'], line_width=4)
p2.line(x='c', y='Q291K', source=source, color=colors['Q291K'], line_width=4)


p.legend.location = 'bottom_right'
p2.legend.location = 'top_left'

# Define the sliders
R_slider = Slider(start=1, end=10000, value=1, step=10, title='Repressors per cell')
M_slider = Slider(start=1, end=2000, value=1, step=10, title='Plasmids per cell')

# Define the javascript callback. 
cb = CustomJS(args={'source':source, 'R_slider':R_slider, 'M_slider':M_slider}, code=js)

# Set up the widgetbox
R_slider.callback = cb
M_slider.callback = cb
box = widgetbox([R_slider, M_slider])


# Define the layout and show plot
layout = bokeh.layouts.gridplot([[p, box],[p2]])
bokeh.io.save(layout)