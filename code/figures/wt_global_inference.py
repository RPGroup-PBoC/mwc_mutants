# -*- coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import mut.stats
import mut.bayes
import mut.thermo
import bokeh.io
import bokeh.plotting
import bokeh.palettes
import bokeh.layouts
import bokeh.transform
from bokeh.models import ColumnDataSource
import imp
imp.reload(mut.bayes)
bokeh.io.output_notebook()
constants = mut.thermo.load_constants()

# Load the data sets
leak_data = pd.read_csv('../../data/csv/Garcia2011_Brewster2014_data.csv')
ind_data = pd.read_csv('../../data/csv/RazoMejia2018_data.csv')
fug_data = pd.read_csv('../../data/csv/fugacity_data.csv')

# Clean up the data and include necessary info
leak_data['IPTGuM'] = 0
leak_data.rename(columns={'repressor':'repressors'}, inplace=True)
ind_data['repressors']  *= 2
ind_data = ind_data[ind_data['repressors'] > 0]
ind_data.rename(columns={'fold_change_A':'fold_change', 'IPTG_uM':'IPTGuM'},
                inplace=True)

# Merge the datasets
ind_data = ind_data[['operator', 'repressors', 'fold_change', 'IPTGuM']]
ind_data['author'] = 'razo-mejia'
leak_data = leak_data[['operator', 'repressors', 'fold_change', 
                       'IPTGuM', 'author']]
data = pd.concat([ind_data, leak_data])
data['idx'] = data.groupby('operator').ngroup() + 1
data.to_csv('../../data/csv/RazoMejia2018_Garcia2011_Brewster2014_tidy.csv', index=False)
#%%
# Load the stan model
model = mut.bayes.StanModel('../stan/global_inference.stan')
# Define the data dictionary.
data_dict = {'N':len(data), 'J':data['idx'].max(), 'idx':data['idx'],
            'fc':data['fold_change'], 'c':data['IPTGuM'], 'R':data['repressors'],
            'ep_ai': 4.5, 'n_ns':4.6E6, 'n':int(2), 'offset':0, 'prefix':0}

# Sample the model
_, samples = model.sample(data_dict)
#%%
# Save the samples to disk along with the summary
summ = model.summarize_parameters(parnames=['ka', 'ki', 'ep_r', 
                                    'ep_r', 'ep_r', 'ep_r'])

# %%
samples['ep_ai'] = 4.5
samples.rename(columns={'ep_r[1]':'ep_r.O1', 'ep_r[2]':'ep_r.O2',
                        'ep_r[3]':'ep_r.O3', 'ep_r[4]':'ep_r.Oid'}, inplace=True)
summ['ep_ai'] = 4.5 
summ.loc[(summ['parameter']=='ep_r') & (summ['dimension']==1), 'operator'] = 'O1'
summ.loc[(summ['parameter']=='ep_r') & (summ['dimension']==2), 'operator'] = 'O2'
summ.loc[(summ['parameter']=='ep_r') & (summ['dimension']==3), 'operator'] = 'O3'
summ.loc[(summ['parameter']=='ep_r') & (summ['dimension']==4), 'operator'] = 'Oid'

samples.to_csv('../../data/csv/razomejia_epai_samples.csv')
summ.to_csv('../../data/csv/razomejia_epai_summary.csv')

# %%
summ
#%% Isolate the posterior medians for each parameter
ka = summ[summ['parameter']=='ka']['median'].values[0]
ki = summ[summ['parameter']=='ki']['median'].values[0]
O1 = summ[(summ['parameter']=='ep_r') & (summ['operator']=='O1')]['median'].values[0]
O2 = summ[(summ['parameter']=='ep_r') & (summ['operator']=='O2')]['median'].values[0]
O3 = summ[(summ['parameter']=='ep_r') & (summ['operator']=='O3')]['median'].values[0]
Oid = summ[(summ['parameter']=='ep_r') & (summ['operator']=='Oid')]['median'].values[0]
ep_ai = np.log(0.7)

# Compute the theoretical curves
c_range = np.logspace(-2, 4, 500)
rep_range = np.logspace(0, 3.5, 500)
ind_dfs = {}
op_dict = {O1:'O1', O2:'O2', O3:'O3'}
for o, v in op_dict.items():
    dfs = {'fold_change':[], 'IPTGuM':[], 'repressors':[]}
    for r in [22, 60, 124, 260, 1740, 1220]:
        fc = mut.thermo.SimpleRepression(R=r, ep_r=o, ka=ka, ki=ki, ep_ai=ep_ai,
                                         effector_conc=c_range).fold_change()
        dfs['fold_change'].append(fc)
        dfs['IPTGuM'].append(c_range)
        dfs['repressors'].append(str(r))

    ind_dfs[v] = dfs

# Isolate the operators
o1_theo = ColumnDataSource(ind_dfs['O1'])
o2_theo = ColumnDataSource(ind_dfs['O2'])
o3_theo = ColumnDataSource(ind_dfs['O3'])

# Compute the Leakiness theory
op_dict[Oid] = 'Oid'
leak_dfs = {'fold_change':[], 'operator':[], 'repressors':[]}
for o, v in op_dict.items():
    fc = mut.thermo.SimpleRepression(R=rep_range, ep_r=o, ka=ka, ki=ki, ep_ai=ep_ai,
                                     effector_conc=0).fold_change()
    leak_dfs['fold_change'].append(fc)
    leak_dfs['repressors'].append(rep_range)
    leak_dfs['operator'].append(v)
leak_theo = ColumnDataSource(leak_dfs)


# Compute the fugacity curves
xs = []
ys = []
cs = []
for n in [64, 52]:
    xs.append(rep_range)
    x = np.exp(-O1)
    pact = 1 / (1 + np.exp(-ep_ai))
    r = pact * rep_range
    B = rep_range * (1 + x) - n * x - 4.6E6
    A = x * (r - n - 4.6E6)
    numer = -B - np.sqrt(B * B - 4 * A * r)
    denom = 2 * A
    ys.append(1 / (1 + (numer / denom) * np.exp(-O1)))
    cs.append(str(n))
fug_theo = ColumnDataSource({'repressors':xs, 'leakiness':ys, 'N':cs})

# %%
# Compute the summarized data
data_sum = data.groupby(['author', 'repressors', 'IPTGuM', 'operator'])['fold_change'].agg(('mean', 'sem')).reset_index()
#%%
# Set up the figure canvases
o1_ax = bokeh.plotting.figure(width=300, height=300, x_axis_label='IPTG [µM]',
                             y_axis_label='fold-change', x_axis_type='log',
                             title=f'Operator O1: {O1:0.1f} kT')
o2_ax = bokeh.plotting.figure(width=300, height=300, x_axis_label='IPTG [µM]',
                             y_axis_label='fold-change', x_axis_type='log',
                             title=f'Operator O2: {O2:0.1f} kT')
o3_ax = bokeh.plotting.figure(width=300, height=300, x_axis_label='IPTG [µM]',
                             y_axis_label='fold-change', x_axis_type='log',
                             title=f'Operator O3: {O3:0.1f} kT')

leak_ax = bokeh.plotting.figure(width=450, height=300, x_axis_label='repressors per cell',
                             y_axis_label='leakiness', x_axis_type='log',
                             y_axis_type='log')
fug_ax = bokeh.plotting.figure(width=450, height=300, x_axis_label='repressors per cell',
                             y_axis_label='leakiness', x_axis_type='log',
                             y_axis_type='log')

# Define the data sources
o1_data = data_sum[(data_sum['operator']=='O1') & (data_sum['author']=='razo-mejia')].copy()
o1_data['repressors'] = o1_data['repressors'].values.astype(int).astype(str)
o2_data = data_sum[(data_sum['operator']=='O2') & (data_sum['author']=='razo-mejia')].copy()
o2_data['repressors'] = o2_data['repressors'].values.astype(int).astype(str)
o3_data = data_sum[(data_sum['operator']=='O3') & (data_sum['author']=='razo-mejia')].copy()
o3_data['repressors'] = o3_data['repressors'].values.astype(int).astype(str)
garcia = data[data['author']=='garcia']
brewster = data[data['author']=='brewster']

# Define the color factors
rep_pal = bokeh.transform.factor_cmap('repressors', palette=bokeh.palettes.Dark2_6, 
                                     factors=['22', '60', '124', '260', '1220', '1740'])
op_pal = bokeh.transform.factor_cmap('operator', palette=bokeh.palettes.Dark2_6, 
                                     factors=['O1', 'O2', 'O3', 'O4'])
n_pal = bokeh.transform.factor_cmap('N', palette=bokeh.palettes.Dark2_6,
    factors=['64', '52'])

# Plot the induction glyphs
o1_ax.circle(x='IPTGuM', y='mean', fill_color='white', legend='repressors',
            source=o1_data, line_color=rep_pal, size=6)
o2_ax.circle(x='IPTGuM', y='mean', fill_color='white', 
            source=o2_data, line_color=rep_pal, size=6)
o3_ax.circle(x='IPTGuM', y='mean', fill_color='white', 
            source=o3_data, line_color=rep_pal, size=6)

# Plot the induction curves
o1_ax.multi_line(xs='IPTGuM', ys='fold_change', color=rep_pal, line_width=1, source=o1_theo)
o2_ax.multi_line(xs='IPTGuM', ys='fold_change', color=rep_pal, line_width=1, source=o2_theo)
o3_ax.multi_line(xs='IPTGuM', ys='fold_change', color=rep_pal, line_width=1, source=o3_theo)

# Plot the leakiness data
leak_ax.circle(x='repressors', y='fold_change', fill_color='white', legend='operator',
              line_color=op_pal, source=garcia, size=6)
leak_ax.diamond(x='repressors', y='fold_change', fill_color='white', 
              line_color=op_pal, source=brewster, size=8)

# plot the leakiness theory
leak_ax.multi_line(xs='repressors', ys='fold_change', color=op_pal, line_width=1, source=leak_theo)

# Plot the fugacity data
fug = fug_data[fug_data['operator']=='O1']
fug['N'] = fug['N'].astype(str)
fug_ax.circle(x='repressor', y='fold_change', color=n_pal, size=6, legend='N',
                source=fug, fill_color='white')

# Plot the fugacity theory
fug_ax.multi_line(xs='repressors', ys='leakiness', line_width=1, source=fug_theo, 
                line_color=n_pal)

# Position the legends
o1_ax.legend.location = 'top_left'
leak_ax.legend.location = 'bottom_left'

# Set up the layout
ind_row = bokeh.layouts.row(o1_ax, o2_ax, o3_ax)
old_row = bokeh.layouts.row(leak_ax, fug_ax)
lay = bokeh.layouts.column(ind_row, old_row)
bokeh.io.show(lay)
#%%
data_sum.head()

#%%
