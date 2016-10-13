import numpy as np
import scipy.special
import scipy.stats as sc
import scipy
import numdifftools as ndt
import datashader as ds
from datashader.bokeh_ext import InteractiveImage
import bokeh.plotting
from fit_bivariate_gaussian_astroML import *

#
# Generic thermodynamic functions
#


# #################
def pact_log(IPTG, ea, ei, epsilon):
    '''
    Returns the probability of a repressor being active as described
    by the MWC model.
    Parameter
    ---------
    IPTG : array-like.
        concentrations of inducer on which to evaluate the function
    ea, ei : float.
        minus log of the dissociation constants of the active and the
        inactive states respectively
    epsilon : float.
        energy difference between the active and the inactive state
    Returns
    -------
    pact : float.
        probability of a repressor of being in the active state.
        Active state is defined as the state that can bind to the DNA.
    '''
    pact = (1 + IPTG * np.exp(ea))**2 / ((1 + IPTG * np.exp(ea))**2 +
                                         np.exp(-epsilon) *
                                         (1 + IPTG * np.exp(ei))**2)

    return pact


# #################
def fold_change_log(IPTG, ea, ei, epsilon, R, epsilon_r):
    '''
    Returns the gene expression fold change according to the
    thermodynamic model with the extension that takes into account the
    effect of the inducer.

    Parameter
    ---------
    IPTG : array-like.
        concentrations of inducer on which to evaluate the function
    ea, ei : float.
        minus log of the dissociation constants of the active and the
        inactive states respectively
    epsilon : float.
        energy difference between the active and the inactive state
    R : array-like.
        repressor copy number for each of the strains. The length of
        this array should be equal to the IPTG array. If only one value
        of the repressor is given it is asssume that all the data points
        should be evaluated with the same repressor copy number
    epsilon_r : array-like
        repressor binding energy. The length of this array
        should be equal to the IPTG array. If only one value of the
        binding energy is given it is asssume that all the data points
        should be evaluated with the same repressor copy number

    Returns
    -------
    fold-change : float.
        gene expression fold change as dictated by the thermodynamic model.
   '''
    return 1 / (1 + 2 * R / 5E6 * pact_log(IPTG, ea, ei, epsilon) *
                (1 + np.exp(-epsilon)) * np.exp(-epsilon_r))


# #################
def fold_change_log_rnap(IPTG, ea, ei, epsilon, R, epsilon_r, P, epsilon_p):
    '''
    Returns the gene expression fold change according to the thermodynamic
    model with the extension that takes into account the effect of the
    inducer.

    Parameter
    ---------
    IPTG : array-like.
        concentrations of inducer on which to evaluate the function
    ea, ei : float.
        minus log of the dissociation constants of the active and the
        inactive states respectively
    epsilon : float.
        energy difference between the active and the inactive state
    R : array-like.
        repressor copy number for each of the strains. The length of this
        array should be equal to the IPTG array. If only one value of the
        repressor is given it is asssume that all the data points should be
        evaluated with
        the same repressor copy number
    epsilon_r : array-like
        repressor binding energy. The length of this array
        should be equal to the IPTG array. If only one value of the binding
        energy is given it is asssume that all the data points
        should be evaluated with the same repressor copy number
    P : float.
        Number of RNAP per cell.
    epsilon_p : float.
        RNAP binding energy
    Returns
    -------
    fold-change : float.
        gene expression fold change as dictated by the thermodynamic model.
   '''
    return (1 + P / 4.6E6 * np.exp(-epsilon_p)) / \
           (1 + P / 4.6E6 * np.exp(-epsilon_p) +
            2 * R / 4.6E6 * pact_log(IPTG, ea, ei, epsilon) *
            (1 + np.exp(-epsilon)) * np.exp(-epsilon_r))


# #################
def bohr_fn(df, ea, ei, epsilon=4.5):
    '''
    Computes the Bohr parameter for the data in a DataFrame df as a
    function of the MWC parameters ea and ei
    Parameters
    ----------
    df : DataFrame
        Pandas DataFrame containing all the data for which to calculate the
        bohr parameter
    ea, ei : float.
        Minus log of the dissociation constants of the active and the
        inactive states respectively.
    epsilon : float.
        energy difference between the active and the inactive state.

    Returns
    -------
    bohr : array-like.
        Array with all the calculated Bohr parameters.
    '''
    bohr_param = []
    for i in range(len(df)):
        pact = pact_log(IPTG=df.iloc[i].IPTG_uM, ea=ea, ei=ei,
                        epsilon=epsilon)
        F = -np.log(2 * df.iloc[i].repressors / 4.6E6 * pact *
                    (1 + np.exp(-4.5)) * np.exp(-df.iloc[i].binding_energy))
        bohr_param.append(F)
    return bohr_param

#
# Non-linear regression
#


# #################
def log_post(param, indep_var, dep_var, epsilon=4.5):
    '''
    Computes the log posterior for a single set of parameters.
    Parameters
    ----------
    param : array-like.
        param[0] = epsilon_a
        ]aram[1] = epsilon_i
    indep_var : n x 3 array.
        series of independent variables to compute the theoretical fold-change.
        1st column : IPTG concentration
        2nd column : repressor copy number
        3rd column : repressor binding energy
    dep_var : array-like
        dependent variable, i.e. experimental fold-change. Then length of this
        array should be the same as the number of rows in indep_var.

    Returns
    -------
    log_post : float.
        the log posterior probability
    '''
    # unpack parameters
    ea, ei = param

    # unpack independent variables
    IPTG, R, epsilon_r = indep_var[:, 0], indep_var[:, 1], indep_var[:, 2]

    # compute the theoretical fold-change
    fc_theory = fold_change_log(IPTG, ea, ei, epsilon, R, epsilon_r)

    # return the log posterior
    return -len(dep_var) / 2 * np.log(np.sum((dep_var - fc_theory)**2))


# #################
def resid(param, indep_var, dep_var, epsilon=4.5):
    '''
    Residuals for the theoretical fold change.

    Parameters
    ----------
    param : array-like.
        param[0] = epsilon_a
        param[1] = epsilon_i
    indep_var : n x 3 array.
        series of independent variables to compute the theoretical
        fold-change.
        1st column : IPTG concentration
        2nd column : repressor copy number
        3rd column : repressor binding energy
    dep_var : array-like
        dependent variable, i.e. experimental fold-change. Then length of
        this array should be the same as the number of rows in indep_var.

    Returns
    -------
    fold-change_exp - fold-change_theory
    '''
    # unpack parameters
    ea, ei = param

    # unpack independent variables
    IPTG, R, epsilon_r = indep_var[:, 0], indep_var[:, 1], indep_var[:, 2]

    # compute the theoretical fold-change
    fc_theory = fold_change_log(IPTG, ea, ei, epsilon, R, epsilon_r)

    # return the log posterior
    return dep_var - fc_theory


# #################
def non_lin_reg_mwc(df, p0,
                    indep_var=['IPTG_uM', 'repressors', 'binding_energy'],
                    dep_var='fold_change_A', epsilon=4.5, diss_const=False):
    '''
    Performs a non-linear regression on the lacI IPTG titration data
    assuming Gaussian errors with constant variance. Returns the parameters
    e_A == -ln(K_A)
    e_I == -ln(K_I)
    and its corresponding error bars by approximating the posterior
    distribution
    as Gaussian.

    Parameters
    ----------
    df : DataFrame.
        DataFrame containing all the titration information. It should at
        minimum contain the IPTG concentration used, the repressor copy
        number for each strain and the binding energy of such strain as the
        independent variables and obviously the gene expression fold-change
        as the dependent variable.
    p0 : array-like (length = 2).
        Initial guess for the parameter values. The first entry is the guess
        for e_A == -ln(K_A) and the second is the initial guess for
        e_I == -ln(K_I).
    indep_var : array-like (length = 3).
        Array of length 3 with the name of the DataFrame columns that contain
        the following parameters:
        1) IPTG concentration
        2) repressor copy number
        3) repressor binding energy to the operator
    dep_var : str.
        Name of the DataFrame column containing the gene expression
        fold-change.
    epsilon : float.
        Value of the allosteric parameter, i.e. the energy difference between
        the active and the inactive state.
    diss_const : bool.
        Indicates if the dissociation constants should be returned instead of
        the e_A and e_I parameters.

    Returns
    -------
    if diss_const  == True:
        e_A : MAP for the e_A parameter.
        de_A : error bar on the e_A parameter
        e_I : MAP for the e_I parameter.
        de_I : error bar on the e_I parameter
    else:
        K_A : MAP for the K_A parameter.
        dK_A : error bar on the K_A parameter
        K_I : MAP for the K_I parameter.
        dK_I : error bar on the K_I parameter
    '''
    df_indep = df[indep_var]
    df_dep = df[dep_var]

    # Extra arguments given as tuple
    args = (df_indep.values, df_dep.values, epsilon)

    # Compute the MAP
    popt, _ = scipy.optimize.leastsq(resid, p0, args=args)

    # Extract the values
    ea, ei = popt

    # Instantiate a numdifftools Hessian object for the log posterior
    hes_fun = ndt.Hessian(log_post)

    # Compute the Hessian at the map
    hes = hes_fun(popt, df_indep.values, df_dep.values)

    # Compute the covariance matrix
    cov = -np.linalg.inv(hes)

    if diss_const:
        # Get the values for the dissociation constants and their
        # respective error bars
        Ka = np.exp(-ea)
        Ki = np.exp(-ei)
        deltaKa = np.sqrt(cov[0, 0]) * Ka
        deltaKi = np.sqrt(cov[1, 1]) * Ki
        return Ka, deltaKa, Ki, deltaKi

    else:
        return ea, cov[0, 0], ei, cov[1, 1]

#
# DataShader scatter plots
#


# #################
def base_plot(df, x_col, y_col, log=False):
    # Define the range to plot chekcing if it is a log scale or not
    if log:
        x_range = (np.min(np.log10(df[x_col])),
                   np.max(np.log10(df[x_col])))
        y_range = (np.min(np.log10(df[y_col])),
                   np.max(np.log10(df[y_col])))
    else:
        x_range = (df[x_col].min(), df[x_col].max())
        y_range = (df[y_col].min(), df[y_col].max())

    # Initialize the Bokeh plot
    p = bokeh.plotting.figure(
        x_range=x_range,
        y_range=y_range,
        tools='save,pan,wheel_zoom,box_zoom,reset',
        plot_width=500,
        plot_height=500)

    # Add all the features to the plot
    p.xgrid.grid_line_color = '#a6a6a6'
    p.ygrid.grid_line_color = '#a6a6a6'
    p.ygrid.grid_line_dash = [6, 4]
    p.xgrid.grid_line_dash = [6, 4]
    p.xaxis.axis_label = x_col
    p.yaxis.axis_label = y_col
    p.xaxis.axis_label_text_font_size = '15pt'
    p.yaxis.axis_label_text_font_size = '15pt'
    p.background_fill_color = '#F4F3F6'
    return p


# #################
def ds_plot(df, x_col, y_col, log=False):
    if log:
        data = np.log10(df[[x_col, y_col]])
    else:
        data = df[[x_col, y_col]]
    p = base_plot(data, x_col, y_col)
    pipeline = ds.Pipeline(data, ds.Point(x_col, y_col))
    return p, pipeline

#
# Automatic gating of the flow cytometry data
#


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

#
# Useful MCMC functions
#


# #################
def hpd(trace, mass_frac):
    """
    Returns highest probability density region given by
    a set of samples.
    CREDIT : Justin Bois BEBi103
    http://bebi103.caltech.edu

    Parameters
    ----------
    trace : array
        1D array of MCMC samples for a single variable
    mass_frac : float with 0 < mass_frac <= 1
        The fraction of the probability to be included in
        the HPD.  For example, `massfrac` = 0.95 gives a
        95% HPD.

    Returns
    -------
    output : array, shape (2,)
        The bounds of the HPD
    """
    # Get sorted list
    d = np.sort(np.copy(trace))

    # Number of total samples taken
    n = len(trace)

    # Get number of samples that should be included in HPD
    n_samples = np.floor(mass_frac * n).astype(int)

    # Get width (in units of data) of all intervals with n_samples samples
    int_width = d[n_samples:] - d[:n-n_samples]

    # Pick out minimal interval
    min_int = np.argmin(int_width)

    # Return interval
    return np.array([d[min_int], d[min_int+n_samples]])

#
# Homoscedastic Gaussian likelihood
#


# #################
def log_likelihood_mcmc(param, indep_var, dep_var, epsilon=4.5):
    """
    Computes the log likelihood probability.
    Parameteres
    -----------
    param : data-frame.
        The parameters to be fit by the MCMC. This must be an array of length
        3 with the following entries
        param[0] = ea == -lnKa
        param[1] = ei == -lnKi
        param[2] = sigma. Homoscedastic error associated with the Gaussian
        likelihood.
    indep_var : n x 3 array.
        series of independent variables to compute the theoretical
        fold-change.
        1st column : IPTG concentration
        2nd column : repressor copy number
        3rd column : repressor binding energy
    dep_var : array-like
        dependent variable, i.e. experimental fold-change. Then length of this
        array should be the same as the number of rows in indep_var.
    epsilon : float.
        Energy difference between the active and inactive state of the
        repressor.

    Returns
    -------
    log_like : float.
        the log likelihood.
    """
    # unpack parameters
    ea, ei, sigma = param

    # unpack independent variables
    IPTG, R, epsilon_r = indep_var.iloc[:, 0], indep_var.iloc[:, 1],\
        indep_var.iloc[:, 2]

    # compute the theoretical fold-change
    fc_theory = fold_change_log(IPTG, ea, ei, epsilon, R, epsilon_r)

    log_like = np.sum((fc_theory - dep_var)**2) / 2 / sigma**2
    return log_like


# #################
def log_post_mcmc(param, indep_var, dep_var, epsilon=4.5,
                  ea_range=[6, -6], ei_range=[6, -6], sigma_range=[0, 1]):
    '''
    Computes the log posterior probability.
    Parameters
    ----------
    param : array-like.
        The parameters to be fit by the MCMC. This must be an array of
        length 3 with the following entries
        param[0] = ea == -lnKa
        param[1] = ei == -lnKi
        param[2] = sigma. Homoscedastic error associated with the Gaussian
        likelihood.
    indep_var : n x 3 array.
        Series of independent variables to compute the theoretical
        fold-change.
        1st column : IPTG concentration
        2nd column : repressor copy number
        3rd column : repressor binding energy
    dep_var : array-like
        Dependent variable, i.e. experimental fold-change. Then length of
        this array should be the same as the number of rows in indep_var.
    ea_range : array-like.
        Range of variables to use in the prior as boundaries for the ea
        parameter.
    ei_range : array-like.
        Range of variables to use in the prior as boundaries for the ei
        parameter.
    sigma_range : array-like.
        Range of variables to use in the prior as boundaries for the sigma
        param.
    epsilon : float.
        Energy difference between the active and inactive state of the
        repressor.
    '''
    # unpack parameters
    ea, ei, sigma = param

    # Set the prior boundaries. Since the variables have a Jeffreys prior, in
    # the log probability they have a uniform prior
    if ea > np.max(ea_range) or ea < np.min(ea_range)\
            or ei > np.max(ei_range) or ei < np.min(ei_range)\
            or sigma > np.max(sigma_range) or sigma < np.min(sigma_range):
        return -np.inf

    return -(len(indep_var) + 1) * np.log(sigma) - \
    log_likelihood_mcmc(param, indep_var, dep_var, epsilon)


# #################
def mcmc_cred_region(IPTG, flatchain, R, epsilon_r,
                     mass_frac=.95, epsilon=4.5):
    '''
    This function takes every element in the MCMC flatchain and computes
    the fold-change for each IPTG concentration returning at the end the
    indicated mass_frac fraction of the fold change.

    Parameters
    ----------
    IPTG : array-like.
        IPTG concentrations on which evaluate the fold change
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
        array of 2 x len(IPTG) with the upper and the lower fold-change HPD
        bound for each IPTG concentration
    '''
    # initialize the array to save the credible region
    cred_region = np.zeros([2, len(IPTG)])

    # loop through IPTG concentrations, compute all the fold changes and
    # save the HPD for each concentration
    for i, c in enumerate(IPTG):
        fc = fold_change_log(c, flatchain[:, 0], flatchain[:, 1], epsilon,
                             R, epsilon_r)
        cred_region[:, i] = hpd(fc, mass_frac)

    return cred_region


# #################
def mcmc_cred_region_rnap(IPTG, flatchain, R, epsilon_r, P, epsilon_p,
                          mass_frac=.95, epsilon=4.5):
    '''
    This function takes every element in the MCMC flatchain and computes
    the fold-change for each IPTG concentration returning at the end the
    indicated mass_frac fraction of the fold change.

    Parameters
    ----------
    IPTG : array-like.
        IPTG concentrations on which evaluate the fold change
    flatchain : array-like.
        MCMC traces for the two MWC parameteres.
        flatchain[:,0] = ea flat-chain
        flatchain[:,1] = ei flat-chain
    R : float.
        Mean repressor copy number.
    epsilon_r : float.
        Repressor binding energy.
    P : float.
        Mean RNAP copy number
    epsilon_p : float.
        RNAP binding energy
    mass_frac : float with 0 < mass_frac <= 1
        The fraction of the probability to be included in
        the HPD.  For example, `massfrac` = 0.95 gives a
        95% HPD.
    epsilon : float.
        Energy difference between active and inactive state.

    Returns
    -------
    cred_region : array-like
        array of 2 x len(IPTG) with the upper and the lower fold-change HPD
        bound for each IPTG concentration
    '''
    # initialize the array to save the credible region
    cred_region = np.zeros([2, len(IPTG)])

    # loop through IPTG concentrations, compute all the fold changes and
    # save the HPD for each concentration
    for i, c in enumerate(IPTG):
        fc = fold_change_log_rnap(c, flatchain[:, 0], flatchain[:, 1], epsilon,
                                  R, epsilon_r, P, epsilon_p)
        cred_region[:, i] = hpd(fc, mass_frac)

    return cred_region


# #################
def mcmc_cred_reg_error_prop(IPTG, flatchain, mass_frac=.95, epsilon=4.5):
    '''
    This function takes every element in the MCMC flatchain and computes
    the fold-change for each IPTG concentration returning at the end the
    indicated mass_frac fraction of the fold change.

    Parameters
    ----------
    IPTG : array-like.
        IPTG concentrations on which evaluate the fold change
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
        array of 2 x len(IPTG) with the upper and the lower fold-change HPD
        bound for each IPTG concentration
    '''
    # initialize the array to save the credible region
    cred_region = np.zeros([2, len(IPTG)])

    # loop through IPTG concentrations, compute all the fold changes and
    # save the HPD for each concentration
    for i, c in enumerate(IPTG):
        fc = fold_change_log(c, flatchain[:, 0], flatchain[:, 1], epsilon,
                             flatchain[:, 2], flatchain[:, 3])
        cred_region[:, i] = hpd(fc, mass_frac)

    return cred_region
