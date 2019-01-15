import bokeh.io
import bokeh.plotting
import bokeh.layouts
import bokeh.palettes
import skimage.io
import skimage.measure
import seaborn as sns
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib
import altair as alt
from scipy.signal import gaussian, convolve

def plotting_style():
    """
    Sets the style to the publication style
    """
    rc = {'axes.facecolor': '#E3DCD0',
          'font.family': 'Lucida Sans Unicode',
          'grid.linestyle': '-',
          'grid.linewidth': 0.5,
          'grid.alpha': 0.75,
          'grid.color': '#ffffff',
          'mathtext.fontset': 'stixsans',
          'mathtext.sf': 'sans',
          'legend.frameon': True,
          'legend.facecolor': '#FFEDCE',
          'figure.dpi': 150}

    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('mathtext', fontset='stixsans', sf='sans')
    sns.set_style('darkgrid', rc=rc)

def personal_style():
    """
    Sets the plotting style to my preference
    """ 
    rc = {'axes.facecolor': '#f1f2f6',  #'#EAECEE',
          'font.family': 'sans-serif',
          'font.style': 'italic',
          'axes.grid': True,
          'axes.edgecolor': 'slategray',
          'axes.spines.right': False,
          'axes.spines.top': False,
          'axes.axisbelow': True,
          'axes.linewidth': 1,
          'axes.titlesize': 9,
          'axes.grid': True,
          'lines.linewidth': 2,
          'lines.dash_capstyle': 'round',
          'grid.linestyle': '-',
          'grid.linewidth': 0.75,
          'grid.alpha': 0.5,
          'grid.color': '#B9BBBC',
          'axes.labelsize': 9,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'legend.fontsize': 8,
          'legend.frameon': False,
          'xtick.color': '#4b4b4b',
          'ytick.color': '#4b4b4b',
          'axes.xmargin': 0.01,
          'axes.ymargin': 0.01,
          'figure.dpi': 100}

    # plt.rc('mathtext', fontset='dejavuserif', sf='sans')
    plt.rc('text.latex', preamble=r'\usepackage{mathpazo}')
    matplotlib.style.use(rc)
    flat = ['#699FCE', '#64767C', '#6B5E86', '#8389B4', '#95B7D8', '#A6DCE8', '#72A2B6',
           '#484B3E', '#6D7960']
    sns.set_palette(flat)
    return flat      
         
def altair_config():
    """Sets the style of Altair plots to the PBoC style"""
    pboc = {'config': {
                'view': {
                     'fill':'#E3DCD1'
                       }, 
                 'axis': {
                     'gridColor':'#ffffff', 
                     'gridWidth':1, 
                     'gridOpacity':0.75,
                     'domain':False, 
                     'tickColor': '#ffffff', 
                     'labelFont':'Lucida Sans Unicode', 
                     'titleFont':'Lucida Sans Unicode',
                     'titleFontWeight': 'normal'
                        },
                 'text':{
                     'font':'Lucida Sans Unicode', 
                     'fontWeight': 'normal'
                         },
                 'legend': {
                     'titleFont': 'Lucida Sans Unicode',
                     'labelFont':'Lucida Sans Unicode', 
                     'titleFontWeight': 'normal'
                       }
                    }
           }
    return pboc

def color_selector(style):
    """
    Select the color palette of your choice.

    Parameters
    ----------
    style: str "mut" or "pboc"
        A string identifier for the style. "mut" gives colors for single and double mutants.
        "pboc" returns the PBoC2e color palette.

    Returns
    -------
    colors: dict
        Dictionary of colors. If "dna", "double", or "inducer" is the selected style,
        keys will be the mutants in upper case. Double mutant keys will be DNA-IND. For
        pboc, the keys will be the typical color descriptors. 

    """
    # Ensure the provided style name makes sense.
    if  style.lower() not in ['mut', 'pboc']:
        raise ValueError("Provided style must be 'pboc' or 'mut'. {} provided.".format(style))

    # Set the color styles and return.
    if  style.lower() == 'mut':
        colors = {'Y20I': '#738FC1', 'Q21A': '#7AA974', 'Q21M': '#AB85AC',
                  'F164T': '#A97C50', 'Q294K': '#5D737E', 'Q294V': '#D56C55',
                  'Q294R': '#B2AF58', 'Y20I-F164T': '#2d98da', 'Y20I-Q294K': '#34495e',
                  'Y20I-Q294V': '#8854d0', 'Q21A-F164T': '#4b6584', 'Q21A-Q294K': '#EE5A24',
                  'Q21A-Q294V': '#009432', 'Q21M-F164T': '#1289A7', 'Q21M-Q294K': '#6F1E51',
                  'Q21M-Q294V': '#006266', 'WT': '#3C3C3C'} 

    elif style.lower() == 'pboc':
        colors = {'green': '#7AA974', 'light_green': '#BFD598',
              'pale_green': '#DCECCB', 'yellow': '#EAC264',
              'light_yellow': '#F3DAA9', 'pale_yellow': '#FFEDCE',
              'blue': '#738FC1', 'light_blue': '#A9BFE3',
              'pale_blue': '#C9D7EE', 'red': '#D56C55', 'light_red': '#E8B19D',
              'pale_red': '#F1D4C9', 'purple': '#AB85AC',
              'light_purple': '#D4C2D9', 'dark_green':'#7E9D90'}
    return colors



def bokeh_traceplot(samples, varnames=None):
    """
    Generate a Bokey traceplot of a series of parameters and their sampling 
    history. As of now, the only shown distribution is 'ecdf'. 

    Parameters
    ----------
    samples: StanFit4Model object
        Sampling output from Stan.
    varnames: list
        List of variable names to include in the plot.
    """
    params = samples.model_pars
    sample_values = samples.extract()
    palette = bokeh.palettes.Category10_10
    
    # Define the standard plot sizes.   
    pairs = []
    if varnames != None:
        iterator = varnames
    else:
        iterator = params
    for p in iterator:
        colors = itertools.cycle(palette)
        if len(np.shape(sample_values[p])) == 1:
            _samples = np.array([sample_values[p]]).T
        else:
            _samples = sample_values[p]
      
        dfs = []
        trace = bokeh.plotting.figure(plot_width=400, plot_height=200, 
                                      x_axis_label='sample number', y_axis_label=f'{p}',
                                     title=f'sampling history for {p}', background_fill_color='#ecedef') 

        dist = bokeh.plotting.figure(plot_width=400, plot_height=200, 
                                     x_axis_label=f'{p}', y_axis_label='ECDF',
                                    title=f'posterior distribution for {p}', background_fill_color='#ecedef')   
        
        # Add visual formatting
        trace.xgrid.grid_line_color = '#FFFFFF'
        trace.ygrid.grid_line_color = '#FFFFFF'
        dist.xgrid.grid_line_color = '#FFFFFF'
        dist.ygrid.grid_line_color = '#FFFFFF'    
        
        for i, color in zip(range(np.shape(_samples)[1]), colors): 
            # Extract the specific values. 
            _values = _samples[:, i]
            x, y = np.sort(_values), np.arange(0, len(_values), 1) / len(_values)
            _df = pd.DataFrame(np.array([x, y, 
                                         _values, np.arange(0, len(_values), 1)]).T, 
                               columns=['x', 'ecdf', 'samples', 'step_no'])            
            dist.line(_df['x'], _df['ecdf'], line_width=2, color=color)
            trace.line(_df['step_no'], _df['samples'], line_width=1, color=color) 
        pairs.append([dist, trace])  
    return bokeh.io.show(bokeh.layouts.gridplot(pairs))