<p align="center">
  <img src="logo.png">
</p>

# The Energetics of Molecular Adaptation in Transcriptional Regulation
[![Build Status](https://travis-ci.org/RPGroup-PBoC/mwc_mutants.svg?branch=publication)](https://travis-ci.org/RPGroup-PBoC/mwc_mutants)

Welcome to the GitHub repository for the MWC mutants project! This repository
contains the entire project history as well as curated scripts to make
reproducibility feasible. 

## Branches
This repository consists of three branches -- `master`, `gh-pages`, and
`publication`. The `master` branch is the primary working branch the authors
used during the research. It contains all project history, processing data
scripts, preliminary figure scripts, and other exploratory data analysis
files. The `gh-pages` brach contains all of the [website
files](https://www.rpgroup.caltech.edu/mwc_mutants) including hosting of the
interactive figures. Finally, the `publication` branch (which you are reading
this on) contains all final processed data, figure generating scripts, data
analysis scripts, and the software module `mut`. If interested in reproducing
the work in this publication, you can execute the scripts on this branch.
Please see the individual directories for more information.

## Installation
To reproduce this work, you will need to use the `mut` module -- a homegrown
Python software package written explicitly for this work. The requirements
can be installed by executing the following command using
[`pip`](pypi.org/project/pip) in the command line:

``` pip install -r requirements.txt ```

The software module itself can be installed locally by executing the command
in the root directory,

``` pip install -e ./ ```


We have written a file `install.sh` that performs both of these steps. To
install, simply run

``` sh install.sh ```

from the root directory. When installed, a new folder `mut.egg-info` will be
installed and is necessary to run any of the code.

## Repository Architecture
This repository is broken up into several directories and subdirectories. Please
see each directory for information regarding each file. 

### **`code`** 
The name says it all. This repository contains all code used in this work to
analyze data and generate figures. It is broken up into several subdirectories which separate the scripts by function.
1. **``analysis`` \|** Contains all Python scripts which perform data
      *analysis* procedures. This includes parameter inference and calibration of various
      statistical models. 
2. **``figures`` \|** Contains all Python scripts which perform data
         *visualization* procedures. This includes generation all main and supplementary
         text figures as well as the interactive figures on the [paper
         website](https://wwww.rpgroup.caltech.edu/mwc_mutants).
3. **``stan`` \|** Contains all inferential models and associated functions
   written in the [Stan probabilistic programming language](http://mc-stan.org).
4. **``processing`` \|** This folder contains a single script
   (`example_processing.py`) that illustrates how processing of the raw data was
   performed. It is kept very generic such that only a few parameters at the top
   of the script need to be changed from experiment to experiment. 

### **`data`**
This directory contains all processed and simulated data used in this work. The
raw flow cytometry data is stored on the [CaltechDATA](http://data.caltech.edu)
research data repository under the [DOI: XXXXXX]().

### **`figures`**
This folder contains `.pdf` files of all main and supplementary text figures. No
code exists in this directory. 

### **`mut`**
This is the heart-and-soul of the repository. It contains a series of Python
files which define the myriad functions used in the processing, analysis, and
visualization of data associated with this work. Pleas navigate to the directory
to see a description of its contents. 

### **`tests`**
We tested a selection of functions we believe are critical to this work. This
includes the gaussian gating routine for processing of flow cytometry,
computation of MWC active probabilities, and scraping of processed datasets to
collect only the accepted experiments. This repository is powered by [Travis
Continuous Integration](http://travis-ci.org) which periodically runs the test
functions to ensure that everything passes. 

## License
![](https://licensebuttons.net/l/by/3.0/88x31.png)

All creative works (writing, figures, etc) are licensed under the [Creative
Commons CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) license. All software is distributed under the standard MIT license as follows

```
Copyright 2019 The Authors 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```