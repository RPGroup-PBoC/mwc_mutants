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
                  'F164T': '#A97C50', 'Q294K': '#EAC264', 'Q294V': '#D56C55',
                  'Q294R': '#919796', 'Y20I-F164T': '#5D737E', 'Y20I-Q294K': '#63AFCC',
                  'Y20I-Q294V': '#9CC7EA', 'Q21A-F164T': '#3C8468', 'Q21A-Q294K': '#63B6AD',
                  'Q21A-Q294V': '#91D0B1', 'Q21M-F164T': '#987D7C', 'Q21M-Q294K': '#A09CB0',
                  'Q21M-Q294V': '#DDBDD4', 'WT': '#3C3C3C'} 
    elif style.lower() == 'pboc':
        colors = {'green': '#7AA974', 'light_green': '#BFD598',
              'pale_green': '#DCECCB', 'yellow': '#EAC264',
              'light_yellow': '#F3DAA9', 'pale_yellow': '#FFEDCE',
              'blue': '#738FC1', 'light_blue': '#A9BFE3',
              'pale_blue': '#C9D7EE', 'red': '#D56C55', 'light_red': '#E8B19D',
              'pale_red': '#F1D4C9', 'purple': '#AB85AC',
              'light_purple': '#D4C2D9', 'dark_green':'#7E9D90'}
    return colors
