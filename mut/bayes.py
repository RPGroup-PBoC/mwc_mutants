# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import pymc3 as pm
import scipy.special
import scipy.optimize
import statsmodels.tools.numdiff as smnd
import theano.tensor as tt


def chains_to_dataframe(fit, var_names=None):
    """
    Converts the generated traces from MCMC sampling to a tidy
    pandas DataFrame.

    Parameters
    ----------
    fit : pystan sampling output
        The raw MCMC output.
    var_names : list of str
        Names of desired parameters. If `None`, all parameters will be
        returned.

    Returns
    -------
    df : pandas DataFrame
        Pandas DataFrame containing all samples from the MCMC.
    """

    data_out = fit.extract()
    if var_names is None:
        _v = var_names = data_out.keys()
        var_names = []
        for v in _v:
            if v != 'lp__':
                var_names.append(v)

    for v in var_names:
        if v not in data_out.keys():
            raise ValueError("Parameter `{}` not found in index.".format(v))

    df = pd.DataFrame([])
    for k in var_names:
        shape = np.shape(data_out[k])
        if len(np.shape(data_out[k])) == 1:
            df.insert(0, k, data_out[k])
        else:
            for n in range(shape[1]):
                df.insert(0, '{}__{}'.format(k, n), data_out[k][:, n])

    # Compute the logp
    logp = [fit.log_prob([data_out[v][i] for v in var_names])
            for i in range(len(df))]
    df.insert(0, 'logp', logp)

    return df


class MarginalizedHomoscedasticNormal(pm.Continuous):
    """
    A bivariate Normal distribution after marginalization of the
    variance, sigma.

    .. math::

        g(\mu, k \vert y) = \left((y - \mu)^2\right)^{-k/2}


    Parameters
    ----------
    mu : PyMC3 RV object
        Mean components of the distribution.
    """

    def __init__(self, mu=None, *args, **kwargs):
        super(MarginalizedHomoscedasticNormal, self).__init__(*args, **kwargs)
        self.mu = mu = pm.theanof.floatX(tt.as_tensor_variable(mu))
        self.median = mu
        self.mode = mu
        self.mean = mu

    def logp(self, values):
        k = values.shape[-1]
        mu = self.mu
        return -0.5 * k * tt.log(tt.sum((values - mu)**2))


class Jeffreys(pm.Continuous):
    """
    Jeffreys prior for a scale parameter.

    Parameters
    ----------
    lower : float, > 0
        Minimum value the variable can take.
    upper : float, > `lower`
        Maximum value the variable can take.
    Returns
    -------
    output : pymc3 distribution
        Distribution for Jeffreys prior.

    Notes
    -----
    This class was adopted from Justin Bois
    github.com/justinbois/bebi103
    """

    def __init__(self, lower=None, upper=None, transform='interval',
                 *args, **kwargs):
        # Check inputs
        if lower is None or upper is None:
            raise RuntimeError('`lower` and `upper` must be provided.')

        if transform == 'interval':
            transform = pm.distributions.transforms.interval(lower, upper)
        super(Jeffreys, self).__init__(transform=transform, *args, **kwargs)
        self.lower = lower = pm.theanof.floatX(tt.as_tensor_variable(lower))
        self.upper = upper = pm.theanof.floatX(tt.as_tensor_variable(upper))

        self.mean = (upper - lower) / tt.log(upper / lower)
        self.median = tt.sqrt(lower * upper)
        self.mode = lower

    def logp(self, value):
        lower = self.lower
        upper = self.upper
        return pm.distributions.dist_math.bound(
            -tt.log(tt.log(upper / lower)) - tt.log(value),
            value >= lower, value <= upper)


def ReparameterizedNormal(name=None, mu=None, sd=None, shape=1):
    """
    A reparameterized (non-centered) normal distribution. This allows for
    more efficient sampling using PyMC3

    Parameters
    ----------
    name :  string
        The name of the RV. The reparameterized version will have this name prepended with "offset_"
    mu : float
        Mean of the normal distribution.
    sd: float
        The standard deviation if the distribtion.
    shape : int
        The shape of the RV. Default is 1

    Returns
    -------
    var : PyMC3 RV object
        The reparameterized distribution.
    """
    if name is None:
        raise RuntimeError("`name` must be provided.")
    if mu is None:
        raise RuntimeError("`mu` must be provided.")
    if sd is None:
        raise RuntimeError("`sd` must be provided.")
    if type(name) is not str:
        raise TypeError(
            "expected type(name) to be string, got {0}.".format(type(name)))

    # Compute the offset.
    offset_var = pm.Normal('offset_{0}'.format(name), mu=0, sd=1, shape=shape)

    # Define the reparameterized variable.
    var = pm.Deterministic(name, mu + offset_var * sd)
    return var
