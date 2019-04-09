---
layout: category
title: Code
permalink: /code/
header:
    overlay_image: /assets/images/seg.png
---

All code used in this work is available through the associated [GitHub
repository](https://www.github.com/rpgroup-pboc/mwc_mutants).


## Jupyter Notebooks  
We have put together several Jupyter notebooks which illustrate the various data
processing and inferential techniques used in the manuscript. 

* [Automated gating of flow cytometry data]()
* [Validation of inferential models]()
* [Inference of the free energy from fold-change data]()

## Stan Models
All inferential models were written using the [Stan probabilistic programming
language](http://mc-stan.org). The models and functions used in this work are
available below:
* [`functions.stan`]()  \| A variety of functions used in this work
* [`DNA_binding_energy.stan`]() \| Stan model for inferring the DNA binding
  energy from a single induction profile.
* [`KaKi_only.stan`]() \| Stan model for inferring the inducer binding constants
  from a single induction profile.
* [`KaKi_epAI.stan]() \| Stan model for inferring the inducer binding constants
  as well as the allosteric energy difference
* [`mean_fold_change.stan]() \| Stan model for inferring the mean fold-change
  and standard deviation from a set of fold-change measurements.

## `mut` Python Module 
