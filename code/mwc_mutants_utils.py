"""
Title:
    mwc_mutants_utils.py
Creation Date:
    2016-07-16
Author(s):
    Manuel Razo-Mejia, Griffin Chure
Purpose:
    This file contains a myriad of functions used for data processing, analysis
    and inference of a variety of data types. These functions are split into
    different categories, which are listed below.

        General Thermodynamic Functions
        -------------------------------
        pact_log:
            Returns the probability of a repressor being active as described
            by the MWC model.

        fold_change_log:
            Returns the gene expression fold change according to the
            thermodynamic model with the extension that takes into account the
            effect of the inducer.

           bohr_fn:
            Computes the Bohr parameter for the data in a DataFrame df as a
            function of the MWC parameters ea and ei.

        Regression and Bayesian Inference
        ---------------------------------
        log_post:
            Computes the log posterior of the fold-change expression for a
            single set of parameters

        resid:
            Residuals for the theoretical fold change.

        non_lin_reg_mwc:
            Performs a non-linear regression on the lacI iptg titration data
            assuming Gaussian errors with constant variance. Returns the
            parameters
                                    e_A = -ln(K_A)
                                    e_I = -ln(K_I)
            and their corresponding error bars by approximating the posterior
            distribution as Gaussian.

        log_likelihood_mcmc:
            Computes the log likelihood probability (homoscedastic gaussian) of
            the fold-change expression for input into the mcmc hammer.

        log_post_mcmc:
            Computes the log posterior probability of the fold-change
            expression for input into the mcmc hammer.

        MCMC Utilities
        --------------
        hpd:
            Returns highest probability density region given by a set of
            samples.

        mcmc_cred_region:
            This function takes every element in the MCMC flatchain and
            computes the fold-change for each iptg concentration returning at
            the end the indicated mass_frac fraction of the fold change.
     
        mcmc_cred_reg_error_prop:
            This function takes every element in the MCMC flatchain and
            computes the fold-change for each iptg concentration returning at
            the end the indicated mass_frac fraction of the fold change.

        Flow Cytometry Data Processing
        ------------------------------
        fit_2D_gaussian:
            This function hacks astroML fit_bivariate_normal to return the
            mean and covariance matrix when fitting a 2D gaussian fuction to
            the data contained in the x_vall and y_val columns of the
            DataFrame df.

        gauss_interval:
            Computes the of the statistic (x - µx)'sum(x - µx) for each of the
            elements in df columns x_val and y_val.

        auto_gauss_gate:
            Function that applies an "unsupervised bivariate Gaussian gate" to
            the data over the channels x_val and y_val.

        Image Processing
        ----------------
        find_zero_crossings:
            This  function computes the gradients in pixel values of an image
            after applying a sobel filter to a given image. This  function is
            later used in the Laplacian of Gaussian cell segmenter
            (log_segmentation) function.

        log_segmentation:
            This function computes the Laplacian of a gaussian filtered image
            and detects object edges as regions which cross zero in the
            derivative.

        average_stack:
            Computes an average image from a provided array of images.

        generate_flatfield:
            Corrects illumination of a given image using a dark image and an
            image of the flat illumination.

        gaussian_subtraction:
            This function applies a gaussian blur to an image and subtracts it
            from the original image. This removes any large-scale
            irregularities in illumination.

        Plotting Configuration
        ----------------------
        set_plotting_style:
            Formats plotting enviroment to that used in Physical Biology of
            the Cell, 2nd edition. To format all plots within a script, simply
            execute `mwc_induction_utils.set_plotting_style() in the preamble.

License: MIT
    Copyright (c) 2017 Rob Phillips group @ California Institute of Technology

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.MIT

"""
import numpy as np
import scipy.special
import scipy.ndimage
import skimage.filters
import skimage.morphology
import skimage.segmentation
import joblib as jlb
import pandas as pd
import scipy.stats as sc
import scipy
import statsmodels.tools.numdiff as smnd # to comput the Hessian matrix
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from fit_bivariate_gaussian_astroML import *

# #################
# Generic thermodynamic functions
# #################

# #################
def pact_log(iptg, ka, ki, epsilon=4.5, n=2):
    '''
    Returns the probability of a repressor being active as described
    by the MWC model.

    Parameter
    ---------
    iptg : array-like.
        Concentrations of inducer on which to evaluate the function.
        All values must be positive.
    ka, ki : float.
        Minus log of the dissociation constants of the active and the
        inactive states respectively.
    epsilon : float.
        Positive log of the energy difference between the active and the
        inactive state.
    n : int
        Number of inducer binding sites.
    Returns
    -------
    pact : float.
        probability of a repressor of being in the active state.
        Active state is defined as the state that can bind to the DNA.

    Raises
    ------
    ValueError
        Thrown if any value of iptg concentration are negative.
    '''
    # Ensure that all values of iptg are positive.
    if (iptg < 0).any():
        raise ValueError('iptg array cannot have negative values.')

    pact = (1 + iptg * np.exp(ka))**n / ((1 + iptg * np.exp(ka))**n +
                                         np.exp(-epsilon) *
                                         (1 + iptg * np.exp(ki))**n)

    return pact


# #################
def fold_change_log(iptg, ka, ki, epsilon, R, epsilon_RA, n=2,
                    quaternary_state=1, nonspec_sites=4.6E6):
    '''
    Returns the gene expression fold change according to the
    thermodynamic model with the extension that takes into account the
    effect of the inducer.

    Parameter
    ---------
    iptg : array-like.
        Concentrations of inducer on which to evaluate the function
    ka, ki : float.
        Minus log of the dissociation constants of the active and the
        inactive states respectively
    epsilon : float.
        Energy difference between the active and the inactive state
    R : array-like.
        Repressor copy number for each of the strains. The length of
        this array should be equal to the iptg array. If only one value
        of the repressor is given it is asssume that all the data points
        should be evaluated with the same repressor copy number
    epsilon_RA : array-like
        Repressor binding energy. The length of this array
        should be equal to the iptg array. If only one value of the
        binding energy is given it is asssume that all the data points
        should be evaluated with the same repressor copy number
    quaternary_state: int
        Prefactor in front of R in fold-change. Default is 2
        indicating that there are two functional heads per repressor molecule.
        This value must not be zero.
    nonspec_sites : int
        Number of nonspecific binding sites in the system.
        This value must be greater than 0.

    Returns
    -------
    fold_change : float.
        Gene expression fold change as dictated by the thermodynamic model.

    Raises
    ------
    ValueError
        Thrown if any entry of the IPTG vector, number of repressors,
        quaternary prefactor, or number of nonspecific binding sites is
        negative. This is also thrown if the quaternary
        state  or number of nonspecific binding sites is 0.


   '''
    # Ensure that IPTG values and R is positive.
    if (iptg < 0).any() or (R < 0).any():
        raise ValueError('iptg and R must be positive.')
    if (quaternary_state <= 0) or (nonspec_sites <= 0):
        raise ValueError('quaternary_state  and nonspec_sites must be greater\
        than zero.')

    return (1 + quaternary_state * R / nonspec_sites *
            pact_log(iptg, ka, ki, epsilon, n) * (1 + np.exp(-epsilon)) *
            np.exp(-epsilon_RA))**-1

# #################
def bohr_fn(df, ka, ki, epsilon=4.5, quaternary_state=2, nonspec_sites=4.6E6):
    '''
    Computes the Bohr parameter for the data in a DataFrame df as a
    function of the MWC parameters ea and ei
    Parameters
    ----------
    df : pandas DataFrame
        Pandas DataFrame containing all the data for which to calculate the
        bohr parameter. The provided dataframe must have a specific structure.
        See notes for more information
    ka, ki : float
        Minus log of the dissociation constants of the active and the
        inactive states respectively.
    epsilon : float
        energy difference between the active and the inactive state.
    quaternary_state: int
        Prefactor in front of R in fold-change. Default is 2
        indicating that there are two functional heads per repressor molecule.
        This value must not be zero.
    nonspec_sites : int
        Number of nonspecific binding sites for the repressor. Default value
        is 4.6E6. This must be larger than zero.

    Returns
    -------
    bohr_param : array-like.
        Array with all the calculated Bohr parameters.

    Raises
    ------
    RuntimeError
        Thrown if provided data frame is empty.
    ValueError
        Thrown if value of nonspecific sites is less than or equal to 0 or if
        quaternary state is negative.

    Notes
    -----
    The provided data frame must have the following column names:

        IPTG Concentration : 'IPTG_uM',
        Number of Repressors : 'repressors'
        Operator Binding Energy : 'binding_energy'

    Note that this function takes the negative log of the binding energy.


    '''

    # Ensure that the length of the data frame is not 0.
    if 0 in np.shape(df):
        raise RuntimeError('df cannot be empty.')
    if quaternary_state < 0:
        raise ValueError('quaternary_state must be positive.')
    if nonspec_sites <= 0:
        raise ValueError('nonspec_sites must be greater than 0.')

    bohr_param = []
    for i in range(len(df)):
        pact = pact_log(iptg=df.iloc[i]['IPTG_uM'], ka=ka, ki=ki,
                        epsilon=epsilon)
        calc_param = -np.log(quaternary_state * df.iloc[i]['repressors'] /
                             nonspec_sites * pact * (1 + np.exp(-epsilon)) *
                             np.exp(-df.iloc[i]['binding_energy']))
        bohr_param.append(calc_param)
    return bohr_param

#============================================================================== 
# Non-linear regression parameter estimation
#============================================================================== 

def nonlin_log_post(param, param_names, indep_var, fc_exp):
    '''
    Computes the log likelihood probability.
    Parameteres
    -----------
    param : array-like.
        The parameters to be fit by the MCMC. For this function there are 4
        possible parameters that the routine can fit:
        1. epsilon_RA
        2. ka == -lnKa
        3. ki == -lnKi
        4. epsilon
        5. R
        Any of these parameters can be fit in any combination as long as the names
        are indicated in in the param_names array. But for these the parameters
        must be fed in the same order as indicated in the param_names array.
    param_names : array-like.
        array containing strings with the names of the parameters given in the
        param array. NOTE: It is important that these names are given in the same
        order as the in the param array.
    indep_var : dictionary.
        series of independent variables to compute the theoretical fold-change.
        These array MUST contain all the variables necessary to compute the
        theoretical fold-change that were not indicated in the param array.
    fc_exp : array-like.
        experimental fold change of each of the data points.

    Returns
    -------
    log_post : float.
        the log posterior probability

    Raises
    ------
    ValueError
        Thrown if quanternary state is negative or if number of nonspecific
        sites is less than or equal to 0.

    RuntimeError
        Thrown if indepvar is not Nx3 or if the length of the dependent
        variables is not the same as the number of elements in the independent
        variables.
    '''

    # Generate dictionary with the given parameters to compute the theoretical
    # fold change.
    fc_param = {**dict(zip(param_names, param)), **indep_var}
    
    # compute the theoretical fold-change
    fc_theory = fold_change_log(**fc_param)

    return - len(indep_var['iptg']) / 2 * np.log(np.sum((fc_exp - fc_theory)**2))

#============================================================================== 

def resid(param, param_names, indep_var, fc_exp):
    '''
    Residuals for the theoretical fold change.

    Returns
    -------
    fold-change_exp - fold-change_theory
    '''
    # Generate dictionary with the given parameters to compute the theoretical
    # fold change.
    fc_param = {**dict(zip(param_names, param)), **indep_var}
    
    # compute the theoretical fold-change
    fc_theory = fold_change_log(**fc_param)

    # return the log posterior
    return fc_exp - fc_theory


# #################
def non_lin_reg(p0, param_names, indep_var, fc_exp):
    '''
    '''
    # Extra arguments given as tuple
    args = (param_names, indep_var, fc_exp)

    # Compute the MAP
    popt, _ = scipy.optimize.leastsq(resid, p0, args=args)

    # Compute the Hessian at the map
    hes = smnd.approx_hess(popt, nonlin_log_post,
                           args=(param_names, indep_var, fc_exp))

    # Compute the covariance matrix
    cov = -np.linalg.inv(hes)

    return popt, cov

# #################
def hpd(trace, mass_frac):
    """
    Returns highest probability density region given by
    a set of samples.
    Parameters
    ----------
    trace : array
        1D array of MCMC samples for a single variable
    mass_frac : float with 0 < mass_frac <= 1
        The fraction of the probability to be included in
        the HPD.  For hreple, `massfrac` = 0.95 gives a
        95% HPD.

    Returns
    -------
    output : array, shape (2,)
        The bounds of the HPD

    Notes
    -----
    We thank Justin Bois (BBE, Caltech) for developing this function.
    http://bebi103.caltech.edu/2015/tutorials/l06_credible_regions.html
    """
    # Get sorted list
    d = np.sort(np.copy(trace))

    # Number of total samples taken
    n = len(trace)

    # Get number of samples that should be included in HPD
    n_samples = np.floor(mass_frac * n).astype(int)

    # Get width (in units of data) of all intervals with n_samples samples
    int_width = d[n_samples:] - d[:n - n_samples]

    # Pick out minimal interval
    min_int = np.argmin(int_width)

    # Return interval
    return np.array([d[min_int], d[min_int + n_samples]])


# #################
# MCMC Utilities
# #################

# #################
def mcmc_cred_region(iptg, flatchain, R, epsilon_r,
                     mass_frac=.95, epsilon=4.5):
    '''
    This function takes every element in the MCMC flatchain and computes
    the fold-change for each iptg concentration returning at the end the
    indicated mass_frac fraction of the fold change.

    Parameters
    ----------
    iptg : array-like.
        iptg concentrations on which evaluate the fold change
    flatchain : array-like.
        MCMC traces for the two MWC parameteres.
        flatchain[:,0] = ea flat-chain
        flatchain[:,1] = ei flat-chain
    R : float.
        Mean repressor copy number.
    epsilon_r : float.
        Repressor binding energy.
    mass_frac : float with 0 < mass_frac <= 1
        The fraction of the probability to be included in
        the HPD.  For example, `massfrac` = 0.95 gives a
        95% HPD.
    epsilon : float.
        Energy difference between active and inactive state.

    Returns
    -------
    cred_region : array-like
        array of 2 x len(iptg) with the upper and the lower fold-change HPD
        bound for each iptg concentration
    '''
    # initialize the array to save the credible region
    cred_region = np.zeros([2, len(iptg)])

    # loop through iptg concentrations, compute all the fold changes and
    # save the HPD for each concentration
    for i, c in enumerate(iptg):
        fc = fold_change_log(c, flatchain[:, 0], flatchain[:, 1], epsilon,
                             R, epsilon_r)
        cred_region[:, i] = hpd(fc, mass_frac)

    return cred_region

# #################
def mcmc_cred_reg_error_prop(iptg, flatchain, mass_frac=.95, epsilon=4.5):
    '''
    This function takes every element in the MCMC flatchain and computes
    the fold-change for each iptg concentration returning at the end the
    indicated mass_frac fraction of the fold change.

    Parameters
    ----------
    iptg : array-like.
        iptg concentrations on which evaluate the fold change
    flatchain : array-like.
        MCMC traces for the two MWC parameteres.
        flatchain[:,0] = ea flat-chain
        flatchain[:,1] = ei flat-chain
        flatchain[:,2] = R flat-chain
        flatchain[:,3] = epsilon_r flat-chain
    mass_frac : float with 0 < mass_frac <= 1
        The fraction of the probability to be included in
        the HPD.  For example, `massfrac` = 0.95 gives a
        95% HPD.
    epsilon : float.
        Energy difference between active and inactive state.

    Returns
    -------
    cred_region : array-like
        array of 2 x len(iptg) with the upper and the lower fold-change HPD
        bound for each iptg concentration
    '''
    # initialize the array to save the credible region
    cred_region = np.zeros([2, len(iptg)])

    # loop through iptg concentrations, compute all the fold changes and
    # save the HPD for each concentration
    for i, c in enumerate(iptg):
        fc = fold_change_log(c, flatchain[:, 0], flatchain[:, 1], epsilon,
                             flatchain[:, 2], flatchain[:, 3])
        cred_region[:, i] = hpd(fc, mass_frac)

    return cred_region

# #################
# Flow Cytometry Data Processing
# #################

# #################
def fit_2D_gaussian(df, x_val='FSC-A', y_val='SSC-A', log=False):
    '''
    This function hacks astroML fit_bivariate_normal to return the mean
    and covariance matrix when fitting a 2D gaussian fuction to the data
    contained in the x_vall and y_val columns of the DataFrame df.

    Parameters
    ----------
    df : DataFrame.
        dataframe containing the data from which to fit the distribution
    x_val, y_val : str.
        name of the dataframe columns to be used in the function
    log : bool.
        indicate if the log of the data should be use for the fit or not

    Returns
    -------
    mu : tuple.
        (x, y) location of the best-fit bivariate normal
    cov : 2 x 2 array
        covariance matrix.
        cov[0, 0] = variance of the x_val column
        cov[1, 1] = variance of the y_val column
        cov[0, 1] = cov[1, 0] = covariance of the data
    '''
    if log:
        x = np.log10(df[x_val])
        y = np.log10(df[y_val])
    else:
        x = df[x_val]
        y = df[y_val]

    # Fit the 2D Gaussian distribution using atroML function
    mu, sigma_1, sigma_2, alpha = fit_bivariate_normal(x, y, robust=True)

    # compute covariance matrix from the standar deviations and the angle
    # that the fit_bivariate_normal function returns
    sigma_xx = ((sigma_1 * np.cos(alpha)) ** 2 +
                (sigma_2 * np.sin(alpha)) ** 2)
    sigma_yy = ((sigma_1 * np.sin(alpha)) ** 2 +
                (sigma_2 * np.cos(alpha)) ** 2)
    sigma_xy = (sigma_1 ** 2 - sigma_2 ** 2) * np.sin(alpha) * np.cos(alpha)

    # put elements of the covariance matrix into an actual matrix
    cov = np.array([[sigma_xx, sigma_xy], [sigma_xy, sigma_yy]])

    return mu, cov

# #################
def gauss_interval(df, mu, cov, x_val='FSC-A', y_val='SSC-A', log=False):
    '''
    Computes the of the statistic
    (x - µx)'sum(x - µx)
    for each of the elements in df columns x_val and y_val.

    Parameters
    ----------
    df : DataFrame.
        dataframe containing the data from which to fit the distribution
    mu : array-like.
        (x, y) location of bivariate normal
    cov : 2 x 2 array
        covariance matrix
    x_val, y_val : str.
        name of the dataframe columns to be used in the function
    log : bool.
        indicate if the log of the data should be use for the fit or not.

    Returns
    -------
    statistic_gauss : array-like.
        array containing the result of the linear algebra operation:
        (x - µx)'sum(x - µx)
    '''

    # Determine that the covariance matrix is not singular
    det = np.linalg.det(cov)
    if det == 0:
        raise NameError("The covariance matrix can't be singular")

    # Compute the vector x defined as [[x - mu_x], [y - mu_y]]
    if log is True:
        x_vect = np.log10(np.array(df[[x_val, y_val]]))
    else:
        x_vect = np.array(df[[x_val, y_val]])

    x_vect[:, 0] = x_vect[:, 0] - mu[0]
    x_vect[:, 1] = x_vect[:, 1] - mu[1]

    # compute the inverse of the covariance matrix
    inv_sigma = np.linalg.inv(cov)

    # compute the operation
    interval_array = np.zeros(len(df))
    for i, x in enumerate(x_vect):
        interval_array[i] = np.dot(np.dot(x, inv_sigma), x.T)

    return interval_array


# #################
def auto_gauss_gate(df, alpha, x_val='FSC-A', y_val='SSC-A', log=False,
                    verbose=False):
    '''
    Function that applies an "unsupervised bivariate Gaussian gate" to the data
    over the channels x_val and y_val.

    Parameters
    ----------
    df : DataFrame.
        dataframe containing the data from which to fit the distribution
    alpha : float. [0, 1]
        fraction of data aimed to keep. Used to compute the chi^2 quantile
        function
    x_val, y_val : str.
        name of the dataframe columns to be used in the function
    log : bool.
        indicate if the log of the data should be use for the fit or not
    verbose : bool.
        indicate if the percentage of data kept should be print

    Returns
    -------
    df_thresh : DataFrame
        Pandas data frame to which the automatic gate was applied.
    '''
    data = df[[x_val, y_val]]
    # Fit the bivariate Gaussian distribution
    mu, cov = fit_2D_gaussian(data, log=log)

    # Compute the statistic for each of the pair of log scattering data
    interval_array = gauss_interval(data, mu, cov, log=log)

    # Find which data points fall inside the interval
    idx = interval_array <= scipy.stats.chi2.ppf(alpha, 2)

    # print the percentage of data kept
    if verbose:
        print('''
        with parameter alpha={0:0.2f}, percentage of data kept = {1:0.2f}
        '''.format(alpha, np.sum(idx) / len(df)))

    return df[idx]

# #################
# Image Processing
# #################

# #################
def ome_split(im):
    """Splits an ome.tiff image into individual channels"""
    if len(np.shape(im)) != 3:
        raise RuntimeError('provided image must be a single image')
    ims = []
    for i in range(np.shape(im)[-1]):
        ims.append(im[:, :, i])
    return ims

# #################
def average_stack(im, median_filt=True):
    """
    Computes an average image from a provided array of images.

    Parameters
    ----------
    im : list or arrays of 2d-arrays
        Stack of images to be filtered.
    median_filt : bool
        If True, each image will be median filtered before averaging.
        Median filtering is performed using a 3x3 square structural element.

    Returns
    -------
    im_avg : 2d-array
        averaged image with a type of int.
    """

    # Determine if the images should be median filtered.
    if median_filt is True:
        selem = skimage.morphology.square(3)
        im_filt = [scipy.ndimage.median_filter(i, footprint=selem) for i in im]
    else:
        im = im_filt

    # Generate and empty image to store the averaged image.
    im_avg = np.zeros_like(im[0]).astype(int)
    for i in im:
        im_avg += i
    im_avg = im_avg / len(im)
    return im_avg


def generate_flatfield(im, im_dark, im_field, median_filt=True):
    """
    Corrects illumination of a given image using a dark image and an image of
    the flat illumination.

    Parameters
    ----------
    im : 2d-array
        Image to be flattened.
    im_dark : 2d-array
        Average image of camera shot noise (no illumination).
    im_field: 2d-array
        Average image of fluorescence illumination.
    median_filt : bool
        If True, the image to be corrected will be median filtered with a
        3x3 square structural element.

    Returns
    -------
    im_flat : 2d-array
        Image corrected for uneven fluorescence illumination. This is performed
        as

        im_flat = ((im - im_dark) / (im_field - im_dark)) *
                   mean(im_field - im_dark)

    Raises
    ------
    RuntimeError
        Thrown if bright image and dark image are approximately equal. This
        will result in a division by zero.
    """

    # Ensure that the same image is not being provided as the bright and dark.
    if np.isclose(im_field, im_dark).all():
        raise RuntimeError('im_bright and im_dark are approximately equal.')

    # Compute the mean difference between the bright and dark image.
    mean_diff = np.mean(im_field - im_dark)

    if median_filt is True:
        selem = skimage.morphology.square(3)
        im_filt = scipy.ndimage.median_filter(im, footprint=selem)
    else:
        im_filt = im

    # Compute and return the flattened image.
    im_flat = ((im_filt - im_dark) / (im_field - im_dark)) * mean_diff
    return im_flat


# #################
def find_zero_crossings(im, selem, thresh):
    """
    This  function computes the gradients in pixel values of an image after
    applying a sobel filter to a given image. This  function is later used in
    the Laplacian of Gaussian cell segmenter (log_segmentation) function. The
    arguments are as follows.

    Parameters
    ----------
    im : 2d-array
        Image to be filtered.
    selem : 2d-array, bool
        Structural element used to compute gradients.
    thresh :  float
        Threshold to define gradients.

    Returns
    -------
    zero_cross : 2d-array
        Image with identified zero-crossings.

    Notes
    -----
    This function as well as `log_segmentation` were written by Justin Bois.
    http://bebi103.caltech.edu/
    """

    # apply a maximum and minimum filter to the image.
    im_max = scipy.ndimage.filters.maximum_filter(im, footprint=selem)
    im_min = scipy.ndimage.filters.minimum_filter(im, footprint=selem)

    # Compute the gradients using a sobel filter.
    im_filt = skimage.filters.sobel(im)

    # Find the zero crossings.
    zero_cross = (((im >= 0) & (im_min < 0)) | ((im <= 0) & (im_max > 0)))\
        & (im_filt >= thresh)

    return zero_cross


# #################
def log_segmentation(im, selem='default', thresh=0.0001, radius=2.0,
                     median_filt=True, clear_border=True, label=False):
    """
    This function computes the Laplacian of a gaussian filtered image and
    detects object edges as regions which cross zero in the derivative.

    Parameters
    ----------
    im :  2d-array
        Image to be processed. Must be a single channel image.
    selem : 2d-array, bool
        Structural element for identifying zero crossings. Default value is
        a 2x2 pixel square.
    radius : float
        Radius for gaussian filter prior to computation of derivatives.
    median_filt : bool
        If True, the input image will be median filtered with a 3x3 structural
        element prior to segmentation.
    selem : 2d-array, bool
        Structural element to be applied for laplacian calculation.
    thresh : float
        Threshold past which
    clear_border : bool
        If True, segmented objects touching the border will be removed.
        Default is True.
    label : bool
        If True, segmented objecs will be labeled. Default is False.

    Returns
    -------
    im_final : 2d-array
        Final segmentation mask. If label==True, the output will be a integer
        labeled image. If label==False, the output will be a bool.

    Notes
    -----
    We thank Justin Bois in his help writing this function.
    https://bebi103.caltech.edu
    """

    # Test that the provided image is only 2-d.
    if len(np.shape(im)) > 2:
        raise ValueError('image must be a single channel!')

    # Determine if the image should be median filtered.
    if median_filt is True:
        selem = skimage.morphology.square(3)
        im_filt = scipy.ndimage.median_filter(im, footprint=selem)
    else:
        im_filt = im
    # Ensure that the provided image is a float.
    if np.max(im) > 1.0:
        im_float = skimage.img_as_float(im_filt)
    else:
        im_float = im_filt

    # Compute the LoG filter of the image.
    im_LoG = scipy.ndimage.filters.gaussian_laplace(im_float, radius)

    # Define the structural element.
    if selem == 'default':
        selem = skimage.morphology.square(3)

    # Using find_zero_crossings, identify the edges of objects.
    edges = find_zero_crossings(im_LoG, selem, thresh)

    # Skeletonize the edges to a line with a single pixel width.
    skel_im = skimage.morphology.skeletonize(edges)

    # Fill the holes to generate binary image.
    im_fill = scipy.ndimage.morphology.binary_fill_holes(skel_im)

    # Remove small objects and objects touching border.
    im_final = skimage.morphology.remove_small_objects(im_fill)
    if clear_border is True:
        im_final = skimage.segmentation.clear_border(im_final, buffer_size=5)

    # Determine if the objects should be labeled.
    if label is True:
        im_final = skimage.measure.label(im_final)

    # Return the labeled image.
    return im_final


# #################
def props_to_df(mask, physical_distance=1, intensity_image=None):
    """
    Converts the output of skimage.measure.regionprops to a nicely
    formatted pandas DataFrame.

    Parameters
    ----------
    mask : 2d-array, int
        Segmentation mask containing objects to be measured.
    physical_distance : int or float
        Interpixel distance of the image. This will be used to
        convert the area measurements to meaningful units.
    intensity_image : 2d-array
        Intensity image for intensity based measurements. If none is
        provided, only region based measurements will be returned.

    Returns
    -------
    df : pandas DataFrame
        Tidy DataFrame containing all measurements.

    """

    # Ensure that there is at least one object in the image.
    if np.max(mask) == 0:
        raise ValueError('no objects found in image.')

    # Define the values that are to be extracted.
    REGIONPROPS = ('area', 'eccentricity', 'solidity',
                   'mean_intensity')

    if intensity_image is None:
        measurements = REGIONPROPS[:-3]
    else:
        measurements = REGIONPROPS

    # Iterate through and extract the props.
    props = skimage.measure.regionprops(mask,
                                        intensity_image=intensity_image)
    for i, p in enumerate(props):
        extracted = []
        for val in measurements:
            extracted.append(p[val])

        if i == 0:
            df = pd.DataFrame(extracted).T
        else:
            df2 = pd.DataFrame(extracted).T
            df = df.append(df2)
    df.columns = measurements
    df['area'] = df['area'] * physical_distance**2
    return df


def example_segmentation(mask, im, bar_length, bounds=True):
    """
    Generates and example segmentation with segmentation mask shown in red over
    the original phase image.

    Parameters
    ----------
    mask : 2d-array, bool
        Boolean mask of segmented objects.
    im : 2d-array, float
        Original image on which the segmentation mask will be overlaid.
    bar_length : int
        Length of scale bar in units of pixels.
    bounds : bool
        If True, only teh bounds of the segmentation mask will be shown around
        each object.

    Returns
    -------
    merge : 3d-array
        Merged segmentation mask.
    """

    # Ensure that the original image is a float and the mask is a bool.
    if np.max(im) > 1:
        im = (im - im.min()) / (im.max() - im.min())
    if np.max(mask) > 0:
        mask = mask > 0

    # Determine if the bounds should be hsown.
    if bounds is True:
        mask = skimage.segmentation.find_boundaries(mask)
    im_copy = np.copy(im)
    im_copy[mask] = 1.0

    return np.dstack((im_copy, im, im))


# #################
# Plotting Configuration
# #################

def ecdf(data):
    """
    Computes the empirical cumulative distribution function (ECDF)
    of a given set of 1D data.

    Parameters
    ----------
    data : 1d-array
        Data from which the ECDF will be computed.

    Returns
    -------
    x, y : 1d-arrays
        The sorted data (x) and the ECDF (y) of the data.
    """

    return np.sort(data), np.arange(len(data))/len(data)


# #################
def set_plotting_style():
    """
    Formats plotting enviroment to that used in Physical Biology of the Cell,
    2nd edition. To format all plots within a script, simply execute
    `mwc_induction_utils.set_plotting_style() in the preamble.
    """
    rc = {'lines.linewidth': 2,
          'axes.labelsize': 18,
          'axes.titlesize': 20,
          'axes.facecolor': '#E3DCD0',
          'xtick.major' : 13,
          'xtick.labelsize': 'large',
          'ytick.labelsize': 13,
          'font.family': 'Lucida Sans Unicode',
          'grid.linestyle': ':',
          'grid.linewidth': 1.5,
          'grid.color': '#ffffff',
          'mathtext.fontset': 'stixsans',
          'mathtext.sf': 'sans',
          'legend.frameon': True,
          'legend.fontsize': 13}
    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('mathtext', fontset='stixsans', sf='sans')
    sns.set_style('darkgrid', rc=rc)
    sns.set_palette("colorblind", color_codes=True)
    sns.set_context('notebook', rc=rc)

def color_dict(names, color_palette='colorblind'):
    '''
    Generates a color palette dictionary given a set of names from the
    indicated seaborn color palette 
    Parameters
    ----------
    names : array-like
        Set of names that will serve as keys for the colors. This can
        be the name of strains or operators or whatever is desired.
    color_palette : str.
        Name of the seaborn palette to be used to generate the dictionary
    '''
    colors = sns.color_palette(color_palette, n_colors=len(names))
    # Change the ugly yellow in the case of the colorblind palette
    if (color_palette=='colorblind') and (len(names) >=5):
        colors[4] = sns.xkcd_palette(['dusty purple'])[0]

    return dict(zip(names, colors))
