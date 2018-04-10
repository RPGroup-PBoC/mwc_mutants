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


def pub_style(return_colors=True):
    """
    Sets the style to the publication style

    Parameters
    ----------
    return_colors: Bool
        If True, this will also return a dictionary

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
    colors = {'green': '#7AA974', 'light_green': '#BFD598',
              'pale_green': '#DCECCB', 'yellow': '#EAC264',
              'light_yellow': '#F3DAA9', 'pale_yellow': '#FFEDCE',
              'blue': '#738FC1', 'light_blue': '#A9BFE3',
              'pale_blue': '#C9D7EE', 'red': '#D56C55', 'light_red': '#E8B19D',
              'pale_red': '#F1D4C9', 'purple': '#AB85AC',
              'light_purple': '#D4C2D9'}
    if return_colors:
        return colors





def _fast_kde(x):
    """
    Note: This code comes directly from the PyMC3 codebase and should only
    be used when specified in the `traceplot` function of this module.

    A fft-based Gaussian kernel density estimate (KDE) for computing
    the KDE on a regular grid.
    The code was adapted from https://github.com/mfouesneau/faststats
    Parameters
    ----------
    x : Numpy array or list
    Returns
    -------
    grid: A gridded 1D KDE of the input points (x).
    xmin: minimum value of x
    xmax: maximum value of x
    """
    x = np.array(x, dtype=float)
    x = x[~np.isnan(x)]
    x = x[~np.isinf(x)]
    n = len(x)
    nx = 200

    # add small jitter in case input values are the same
    x += np.random.uniform(-1E-12, 1E-12, size=n)
    xmin, xmax = np.min(x), np.max(x)

    # compute histogram
    bins = np.linspace(xmin, xmax, nx)
    xyi = np.digitize(x, bins)
    dx = (xmax - xmin) / (nx - 1)
    grid = np.histogram(x, bins=nx)[0]

    # Scaling factor for bandwidth
    scotts_factor = n ** (-0.2)
    # Determine the bandwidth using Scott's rule
    std_x = np.std(xyi)
    kern_nx = int(np.round(scotts_factor * 2 * np.pi * std_x))

    # Evaluate the gaussian function on the kernel grid
    kernel = np.reshape(gaussian(kern_nx, scotts_factor * std_x), kern_nx)

    # Compute the KDE
    # use symmetric padding to correct for data boundaries in the kde
    npad = np.min((nx, 2 * kern_nx))

    grid = np.concatenate([grid[npad: 0: -1], grid, grid[nx: nx - npad: -1]])
    grid = convolve(grid, kernel, mode='same')[npad: npad + nx]

    norm_factor = n * dx * (2 * np.pi * std_x ** 2 * scotts_factor ** 2) ** 0.5

    grid = grid / norm_factor

    return grid, xmin, xmax