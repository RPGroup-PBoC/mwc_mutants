# `mut`

Welcome to the computational guts of the paper. This module contains all
functions used for the processing, analysis, and visualization of the data in
this work. All functions are documented. The following is a brief description of
each module. 

* `__init__.py` \| Standard initiation file for the module. 
* `_fit_bivariate_normal_AstroML.py` \| A copied section from the
  [AstroPy](https://www.astropy.org/) Python package. We used one of the
  functions to fit a bivariate normal distribution to flow cytometry data for an
  automated gating routine. 
* `bayes.py` \| A collection of functions for performing Bayesian inference via
  Stan using the [PyStan](https://pystan.readthedocs.io/en/latest/) interface.
  This file  contains a class `StanModel`.
* `flow.py` \| Contains the automatic gating function and associated helper
  functions from AstroPy.
  which streamlines a series of PyStan MCMC sampling operations. 
* `io.py` \| Input-output functions for scraping through individual experiments
  and reading the metadata to determine if the experiment is accepted or
  rejected. 
* `stats.py` \| Various functions for computing statistics from Pandas
  DataFrames. 
* `thermo.py` \| Functions related to the thermodynamic model of the inducible
  simple repression motif. It is composed of two classes, `MWC` and
  `SimpleRepression`. A function `load_constants` returns a dictionary
  containing the values of the parameters for the wild-type repressor. 
* `viz.py` \| Functions and color definitions for the plotting style. 