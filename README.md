# Quantifying Uncertainty in Subsystem Systems

This repository contains the companion code repository for [Quantifying Uncertainty in Subsurface Systems](AmazonLink) by Céline Scheidt, Lewis Li, and Jef Caers (John Wiley & Sons). 

### About The Book
Summary of book...

### About This Repository
This repository implements the various UQ strategies discussed in the book. The source code for the algorithms can be found under the [src](lewisli/QUSS/src/) folder. The codebase was developed in [MATLAB](https://www.mathworks.com/)

For illustrative purposes, a set of [Jupyter](http://jupyter.org/) tutorials have been prepared. They are as follows:

1. [Dimension Reduction](tutorials/DimRed.ipynb): Showcase of various dimension reduction techniques discussed in Chapter 3.
2. [DGSA](ltutorials/DGSA.ipynb): Implementation of Distance Based Sensitivty Analysis from Chapter 4.
3. [Bayesian Evidential Learning](tutorials/DGSA.ipynb). Methodology discussed in Chapter 7, implemented for the Libyan Oil Reservoir case.
4. [SIR](tutorials/SIR.ipynb): The Sequential Importance Resampling methodology from Chapter 7 applied to the same Libyan Oil Reservoir.

These tutorials can be viewed directly in the browser, or download and re-run. 

#### Installation

To re-run the examples, the following dependencies must be met:

* [MATLAB](https://www.mathworks.com/)
* [MATLAB Engine API for Python](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html)
* [jupyter](http://jupyter.org/) or can be installed as part of [Anaconda](http://softwaremaniacs.org/soft/highlight/en/) package
* [matlab_kernel](https://github.com/Calysto/matlab_kernel) for interfacing MATLAB with jupyter

Once installed, jupyter can be started from the command line with

	jupyter notebook
    
Navigate to the *tutorials* folder and select appropiate tutorial to load.

### Licensing
This repository is released under the MIT License.



