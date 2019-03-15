# -*- coding: utf-8 -*-    
import sys
import pandas as pd
from chrysalids.model import mutant_app
from chrysalids.wildtype import induction_app
from chrysalids.dna_muts import dna_app
import bokeh.io
import bokeh.plotting
from bokeh.models import ColumnDataSource, Div
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select, Slider, RadioButtonGroup
import bokeh.themes
import glob
import  mut.viz
colors = mut.viz.color_selector('pboc')

# ###########################
# DATA LOADING
# ############################
ind_data = pd.read_csv('mwc_mutants/data/csv/RazoMejia2018_data.csv')
ind_data['repressors'] *= 2
ind_data.rename(columns={'fold_change_A':'fold_change',
                        'IPTG_uM':'IPTGuM'}, inplace=True)
ind_data = ind_data[ind_data['repressors'] > 0 ].copy()

# ###########################
# MUTANT DATA PRUNING
# ##########################
mut_data = pd.read_csv('mwc_mutants/data/csv/compiled_data.csv')
dbohr_stats = pd.read_csv('mwc_mutants/data/csv/empirical_F_statistics.csv')
epRA_stats = pd.read_csv('mwc_mutants/data/csv/DNA_binding_energy_summary.csv')
DNA_data = mut_data[((mut_data['class']=='WT') | (mut_data['class']=='DNA')) & 
                    (mut_data['operator']=='O2')]
DNA_stats = dbohr_stats[dbohr_stats['class']=='DNA'].copy()

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

print(epRA_stats.head())
theme = bokeh.themes.Theme(json=theme_json)
tab1 = mutant_app()
tab2 = induction_app(ind_data)
tab3 = dna_app(DNA_data, DNA_stats, epRA_stats)
tabs = bokeh.models.widgets.Tabs(tabs=[tab1, tab2, tab3])
bokeh.io.curdoc().theme = theme
bokeh.io.curdoc().add_root(tabs)