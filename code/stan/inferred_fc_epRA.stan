/* 
* Inferring Fold-Change and DNA Binding Energy
* ------------------------------------------------------------------------------
* Author: Griffin Chure
* License: MIT
* 
* Description
* ------------------------------------------------------------------------------
* This model samples the posterior distribution for the DNA binding energy from 
* a single induction profile of a DNA binding mutant. This model assumes that
* the fold-change measurements are drawn from a normal distribution with a mean 
* mu and a standard deviation sigma. The DNA binding energy is also assumed to
* be drawn from a normal distribution. In this model, the inferred mean of the 
* the fold-change in gene expression is restricted to the mean