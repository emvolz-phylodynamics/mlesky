
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mlesky : Maximum likelihood inference of effective population size through time using GMRF-skygrid approach

This package is related to previous Bayesian implementations of the
Bayesian skygrid model with the following notable differences:

  - The GMRF process takes place on the 2nd-order difference of log(Ne),
    which is more similar to the `skygrowth` model of Volz & Didelot
  - When computing the GMRF likelihood, the smoothing parameter (ie the
    precision of random walk) is fixed
  - We use a novel cross-validation approach for selecting the smoothing
    parameter

## Installation

In R, install the `devtools` package and run

    devtools::install_github('emvolz-phylodynamics/mleksy')

## Documentation

See the vignettes in the R package for examples of how to use mlesky.
