import sys
sys.path.insert(0, '../../../')
import numpy as np
import pandas as pd
from bokeh.themes import Theme
import mut.viz
import mut.thermo
import bokeh.io
import bokeh.plotting
from bokeh import events
from bokeh.models import ColumnDataSource, Div, CustomJS, CDSView, IndexFilter
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select, Slider, RadioButtonGroup, Button,  RangeSlider
from bokeh.transform import factor_cmap, factor_mark

pboc = mut.viz.color_selector('pboc')
constants  = mut.thermo.load_constants()
bokeh.plotting.output_file("DNA_mutants.html")

# Load the necessary data sets and restrict to DNA mutants 
fc_data = pd.read_csv('../../../data/csv/compiled_data.csv')
epRA_data = pd.read_csv('../../../data/csv/DNA_binding_energy_summary.csv')
DNA_data = fc_data[fc_data['class']=='DNA'].copy()

# Define constants for theory plotting
c_range = np.logspace(-2, 4, 200)
bohr_range = np.linspace(-8, 8, 200)

#  Compute the bohr using the different DNA binding energies
for r in DNA_data['repressors'].unique():
    for g, d in DNA_data.groupby(['mutant', 'repressors']):
        param = epRA_data[(epRA_data['parameter']=='ep_RA') & 
                         (epRA_data['mutant']==g[0]) & 
                         (epRA_data['repressors']==r)][['median', 'hpd_min', 
                                                        'hpd_max']]
        _c_data, _ep = np.meshgrid(d['IPTGuM'], param)
        arch = mut.thermo.SimpleRepression(R=g[1], ep_r=_ep, ka=139, ki=0.53,
                                           ep_ai=4.5, n_sites=2, 
                                           effector_conc=_c_data).bohr_parameter()
        DNA_data.loc[(DNA_data['mutant']==g[0]) & 
                     (DNA_data['repressors']==g[1]), f'bohr_R{int(r)}_fit'] = arch[0, :]

# Initialize the fit strain identifier
DNA_data.loc[DNA_data['repressors'] == 260, 'fit_strain'] = 1
DNA_data.loc[DNA_data['repressors'] != 260, 'fit_strain'] = 0

# Convert repressor ID to string for color mapping
DNA_data['repressors'] = DNA_data['repressors'].astype('str')

# Define the source data sets. 
Y20I_source = ColumnDataSource(DNA_data[DNA_data['mutant']=='Y20I'].groupby(
                ['repressors', 'IPTGuM', 'fit_strain']).mean().reset_index())
Q21A_source = ColumnDataSource(DNA_data[DNA_data['mutant']=='Q21A'].groupby(
                ['repressors', 'IPTGuM']).mean().reset_index())
Q21M_source = ColumnDataSource(DNA_data[DNA_data['mutant']=='Q21M'].groupby(
                ['repressors', 'IPTGuM']).mean().reset_index())


# ##############################################################################
# INDUCTION FIGURE AXES
# ##############################################################################
Y20I_ind = bokeh.plotting.figure(width=300, height=300, x_axis_type='log',
                                x_axis_label='IPTG [µM]', 
                                y_axis_label='fold-change',
                                title='Y20I')
Q21A_ind = bokeh.plotting.figure(width=300, height=300, x_axis_type='log',
                                x_axis_label='IPTG [µM]', 
                                y_axis_label='fold-change',
                                title='Q21A') 
Q21M_ind = bokeh.plotting.figure(width=300, height=300, x_axis_type='log',
                                x_axis_label='IPTG [µM]', 
                                y_axis_label='fold-change',
                                title='Q21M') 
# ##############################################################################
# COLLAPSE FIGURE AXES
# ##############################################################################
Y20I_collapse = bokeh.plotting.figure(width=300, height=300,
                                      x_axis_label='free energy [kT]',
                                      y_axis_label='fold_change',
                                      title='Y20I')
Q21A_collapse = bokeh.plotting.figure(width=300, height=300,
                                      x_axis_label='free energy [kT]',
                                      y_axis_label='fold_change',
                                      title='Y20I')
Q21M_collapse = bokeh.plotting.figure(width=300, height=300,
                                      x_axis_label='free energy [kT]',
                                      y_axis_label='fold_change',
                                      title='Q21M')
# ##############################################################################
# GLYPH POPULATION
# ##############################################################################
REPS = list(DNA_data['repressors'].unique().astype('str'))
FIT = ['0', '1']
MARK = ['circle', 'circle_x']
Y20I_ind.scatter(x='IPTGuM', y='fold_change', source=Y20I_source,
                color=factor_cmap('repressors', 'Category10_4', REPS),
                marker=factor_mark('fit_strain', MARK, FIT))
Y20I_collapse.circle(x='bohr_R260_fit', y='fold_change', source=Y20I_source,
                color=factor_cmap('repressors', 'Category10_4', REPS))
Q21A_ind.circle(x='IPTGuM', y='fold_change', source=Q21A_source,
            color=factor_cmap('repressors', 'Category10_4', REPS))
Q21A_collapse.circle(x='bohr_R260_fit', y='fold_change', source=Q21A_source,
            color=factor_cmap('repressors', 'Category10_4', REPS))
Q21M_ind.circle(x='IPTGuM', y='fold_change', source=Q21M_source,
            color=factor_cmap('repressors', 'Category10_4', REPS))
Q21M_collapse.circle(x='bohr_R260_fit', y='fold_change', source=Q21M_source,
           color=factor_cmap('repressors', 'Category10_4', REPS))

Y20I_source.data
# ##############################################################################
# STATIC THEORY CURVES
# ##############################################################################
ref_arch = mut.thermo.SimpleRepression(R=260, ep_r=-13.9, ka=139, ki=0.53, 
                                       n_sites=2, ep_ai=4.5, 
                                       effector_conc=c_range).fold_change()
Y20I_ind.line(c_range, ref_arch, color='gray', line_dash='dotted', line_width=2)
Q21A_ind.line(c_range, ref_arch, color='gray', line_dash='dotted', line_width=2)
Q21M_ind.line(c_range, ref_arch, color='gray', line_dash='dotted', line_width=2)

collapse = (1 + np.exp(-bohr_range))**-1
Y20I_collapse.line(bohr_range, collapse, color='black', line_width=2)
Q21A_collapse.line(bohr_range, collapse, color='black', line_width=2)
Q21M_collapse.line(bohr_range, collapse, color='black', line_width=2)

epRA_data[(epRA_data['repressors']==260) & (epRA_data['parameter']=='ep_RA')]
# ##############################################################################
# FORMATTING AND LAYOUT
# ##############################################################################
layout = bokeh.layouts.gridplot([[Y20I_ind, Q21A_ind, Q21M_ind],
                                 [Y20I_collapse, Q21A_collapse, Q21M_collapse]])

Y20I_source.data
# ##############################################################################
# THEME INFORMATION
# ##############################################################################
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
bokeh.io.show(layout)
