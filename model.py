import sys
sys.path.insert(0, '../../')
import mut.viz
import mut.thermo
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
import itertools
from bokeh.models import ColumnDataSource, Div, CustomJS
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select, Slider, RadioButtonGroup, Button
import glob
constants = mut.thermo.load_constants()
pboc = mut.viz.color_selector('pboc')

c_range = np.logspace(-3, 3, 500)

ref_source = ColumnDataSource(data=dict(c=c_range, fc=c_range))
mut_source = ColumnDataSource(data=dict(c=c_range, fc=c_range))
mut_points = ColumnDataSource(data=dict(c=c_range, fc=c_range))
bokeh.plotting.output_file("model.html")

# ref_fc_callback.args['epRA'] = ref_epRA_slider
p = bokeh.plotting.figure(width=400, height=400, x_axis_type='log')
p.line(x='c', y='fc', source=ref_source, color='black')

ref_callback  = CustomJS(args=(dict(source=ref_source)), code="""
                var data = source.data;
                var ep = cb_obj.value;
                var c = data['c'];
                var fc = data['fc'];

                for (var i = 0; i< c.length; i++) {
                                    fc[i] = ep * c[i]; 
                }
                source.change.emit()
                """)
ref_epRA_slider = Slider(start=-20, end=-2, value=-10, step=0.1,
title='Repressors', callback=ref_callback)
layout = bokeh.layouts.row(ref_epRA_slider, p)
bokeh.io.save(layout)
# var pact = (Math.pow(1 + c[i]/139, 2) / (Math.pow(1 + c[i]/139, 2) + Math.exp(-4.5) * Math.pow(1 + c[i]/0.53, 2));
