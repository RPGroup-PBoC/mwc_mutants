# -*- coding: utf-8 -*-    
import sys
from chrysalids.model import mutant_app
import bokeh.io
import bokeh.plotting
from bokeh.models import ColumnDataSource, Div
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select, Slider, RadioButtonGroup
import bokeh.themes
import glob
import  mut.viz
colors = mut.viz.color_selector('pboc')
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
            'minor_tick_line_color': "white",
            'axis_label_text_font': 'Lucida Sans',
            'axis_label_text_font_style': 'normal'
            },
            'Grid': {
                'grid_line_color': None,
            },
            'Text': {
                'text_font_style': 'normal',
               'text_font': 'Lucida Sans'
            },
            'Title': {
                'background_fill_color': '#FFEDC0',
                'text_font_style': 'normal',
                'align': 'center',
                'text_font': 'Lucida Sans',
                'offset': 2,
            }}}
theme = bokeh.themes.Theme(json=theme_json)
tab1 = mutant_app()
tabs = bokeh.models.widgets.Tabs(tabs=[tab1])
bokeh.io.curdoc().theme = theme
bokeh.io.curdoc().add_root(tabs)