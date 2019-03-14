# -*- coding: utf-8 -*-    
import sys
from chrysalids.chrysalids import mutant_app
import bokeh.io
import bokeh.plotting
from bokeh.models import ColumnDataSource, Div
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select, Slider, RadioButtonGroup
import glob

tab1 = mutant_app()
tabs = bokeh.models.widgets.Tabs(tabs=[tab1])
bokeh.io.curdoc().add_root(tabs)