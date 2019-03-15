import sys
sys.path.insert(0, '../../')
import mut.viz
import mut.thermo
import numpy as np
import pandas as pd
import bokeh.io
import bokeh.plotting
import itertools
from bokeh.models import ColumnDataSource, Div
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select, Slider, RadioButtonGroup, Button, Dropdown
import glob
import os
constants = mut.thermo.load_constants()
pboc = mut.viz.color_selector('pboc')

def dna_app(data, stats, fit_stats):
    data = data.copy()
    summary = data.groupby(['class',
                                 'mutant', 
                                 'repressors', 
                                 'IPTGuM']).agg(('mean', 'sem')).reset_index()
    # Define parameter ranges
    c_range = np.logspace(-2, 4, 500)
    bohr_range = np.linspace(-20, 20, 500)

    # Compute the reference wild-type profile. 
    ref_fc = mut.thermo.SimpleRepression(R=260, ep_r=-13.9,
                                        ka=139, ki=0.53,
                                        ep_ai=4.5, effector_conc=c_range).fold_change()

    # Compute the reference bohr parameter. 
    ref_bohr = mut.thermo.SimpleRepression(R=260, ep_r=-13.9,
                                        ka=139, ki=0.53,
                                        ep_ai=4.5, effector_conc=data['IPTGuM']).bohr_parameter()
    data['ref_bohr'] = ref_bohr
    
    # #######################
    # INSTANTIATING FIGURES 
    # ####################### 
    p_Y20I = bokeh.plotting.figure(width=400, height=300, 
                                x_axis_type='log',
                                y_axis_label='fold-change',
                                x_axis_label='IPTG [µM]',
                                title='    MUTANT Y20I    ',
                                y_range=[-0.1, 1.1])
    p_Y20I_delF = bokeh.plotting.figure(width=400, height=300, 
                                x_axis_type='log',
                                y_axis_label='∆F [kT]',
                                x_axis_label='IPTG [µM]',
                                y_range=[-8, 8],
                                x_range=[1E-2, 1E4])
    p_Q21A = bokeh.plotting.figure(width=400, height=300, 
                                x_axis_type='log',
                                y_axis_label='fold-change',
                                x_axis_label='IPTG [µM]',
                                title='    MUTANT Q21A   ',
                                y_range=[-0.1, 1.1])
    p_Q21A_delF = bokeh.plotting.figure(width=400, height=300, 
                                x_axis_type='log',
                                y_axis_label='∆F [kT]',
                                x_axis_label='IPTG [µM]',
                                y_range=[-8, 8],
                                x_range=[1E-2, 1E4])
    p_Q21M = bokeh.plotting.figure(width=400, height=300, 
                                x_axis_type='log',
                                y_axis_label='fold-change',
                                x_axis_label='IPTG [µM]',
                                title='   MUTANT Q21M    ',
                                y_range=[-0.1, 1.1])
    p_Q21M_delF = bokeh.plotting.figure(width=400, height=300, 
                                x_axis_type='log',
                                y_axis_label='∆F [kT]',
                                x_axis_label='IPTG [µM]',
                                y_range=[-8, 8],
                                x_range=[1E-2, 1E4])
    # Define the selector for repressor copy number
    fit_selector = Select(title='Fit Strain', value='R = 260',
                         options=['R = 60 / cell',
                                  'R = 124 / cell',
                                  'R = 260 / cell',
                                  'R = 1220 / cell'])

    # ##########################
    # SOURCE DATA DEFINITIONS 
    # ##########################
    # There is surely a better way to do this, but this is the most obvious
    Y20I_source = ColumnDataSource(dict(c=[], R60=[], R124=[], R260=[], R1220=[]))
    Q21A_source = ColumnDataSource(dict(c=[], R60=[], R124=[], R260=[], R1220=[]))
    Q21M_source = ColumnDataSource(dict(c=[], R60=[], R124=[], R260=[], R1220=[]))
    source_dict = {'Y20I':Y20I_source, 'Q21A':Q21A_source, 'Q21M':Q21M_source}

    # ##########################
    # DATA PLOTS 
    # ##########################
    axes = {'Y20I':[p_Y20I, p_Y20I_delF],
            'Q21A':[p_Q21A, p_Q21A_delF],
            'Q21M':[p_Q21M, p_Q21M_delF]}
    rep_colors = {60:pboc['blue'], 124:pboc['green'],
                  260:pboc['purple'], 1220:pboc['red']}

    def _update_fit_source():
        # Determine the epRA medians
        R = fit_selector.value
        R = int(R.split('=')[-1].split('/')[0])
        _fit_stats = fit_stats[(fit_stats['repressors']==R) &\
                               (fit_stats['operator']=='O2') &\
                               (fit_stats['parameter']=='ep_RA')]
        for g, d in _fit_stats.groupby('mutant'):
            epRA_median = d['median'].values[0]
            _dict = {'c':c_range}
            for r in data['repressors'].unique():
                arch = mut.thermo.SimpleRepression(R=r, ep_r=epRA_median,
                                                   ka=139, ki=0.53, ep_ai=4.5,
                                                   n_sites=2, 
                                                   effector_conc=c_range).fold_change()
                _dict[f'R{int(r)}'] = arch
            source_dict[g].data = _dict

    def _plot_foldchange_data(data, group, **kwargs):
        # Assign axis
        fc_ax = axes[group[0]][0]

        # Plot mutant fold-change data
        fc_ax.circle(data['IPTGuM'], data['fold_change']['mean'], 
                        color=rep_colors[group[1]])

        # Compute errorbars
        xerrs = [(data.iloc[i]['IPTGuM'], 
                  data.iloc[i]['IPTGuM']) for i in range(len(data))]
        yerrs = [(data.iloc[i]['fold_change']['mean']-\
                  data.iloc[i]['fold_change']['sem'], 
                  data.iloc[i]['fold_change']['mean']+\
                  data.iloc[i]['fold_change']['sem']) for i in range(len(data))]
        fc_ax.multi_line(xerrs, yerrs, line_width=1, color=rep_colors[g[1]])

    def _plot_foldchange_fit(group):
        fc_ax = axes[group[0]][0]
        _update_fit_source()
        # Plot WT theory curve
        fc_ax.line(x='c', y=f'R{int(group[1])}', source=source_dict[group[0]], line_width=2, 
                    color=rep_colors[group[1]])

    def _plot_deltaF(data, summary='median'):
        dbohr  = data[data['parameter']=='delta_bohr_corrected']
        _ax = axes[g[0]][1]
        _ax.circle(dbohr['IPTGuM'], dbohr[summary], color=rep_colors[g[1]],
                    size=5)
        _ax.line(dbohr['IPTGuM'], dbohr[summary], color=rep_colors[g[1]], line_width=1)

        xs = [(x, x) for x in dbohr['IPTGuM'].values]
        ys = [ (dbohr.iloc[i]['hpd_min'], 
                dbohr.iloc[i]['hpd_max']) for i in range(len(dbohr))]
        _ax.multi_line(xs, ys, line_width=1.5,
        color=rep_colors[g[1]])

    def _iterate_foldchange():
        for g, _ in data[data['class']=='DNA'].groupby(['mutant', 'repressors']):
            _plot_foldchange_fit(fit_stats, g)
    
    # ##########################
    # DELTA F DATA
    # ##########################
    for g, d in summary[summary['class']=='DNA'].groupby(['mutant', 'repressors']):
        _plot_foldchange_data(d, g)
        _plot_foldchange_fit(g)

    for g, d in stats.groupby(['mutant', 'repressors']):
        _plot_deltaF(d)

    for a, b in axes.items():
        b[1].line(c_range, np.zeros_like(c_range), color='black', 
                   line_width=2)

    # Configure selectors
    fit_selector.on_change('value', lambda attr, old, new: _update_fit_source())
    lay = layout([[fit_selector], [p_Y20I, p_Q21A, p_Q21M],
                   [p_Y20I_delF,  p_Q21A_delF, p_Q21M_delF]])
    tab = bokeh.models.widgets.Panel(child=lay, title='DNA Binding Mutants')
    return tab