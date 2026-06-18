# multipointR

[![R-CMD-check](https://github.com/mjemons/multipointR/actions/workflows/R-CMD-check.yml/badge.svg?branch=main)](https://github.com/mjemons/multipointR/actions/workflows/R-CMD-check.yml)

`multipointR` is a package to compare the distribution of cells in an image
or across multiple images with labelled point process models.

On a single image, a point process model (`ppm`) 
can be used to model the spatial distribution of a cell type point pattern as a function of spatial covariates
while accounting for natural spacing of cells. The main model class considered in `multipointR` 
are inhomogeneous Gibb's point process models.

Across multiple images, users can either compare coefficients from separately fitting multiple univariate `ppm` models, or 
fit one joint model across all images with `mppm`. 

`multipointR` provides an interface between `SpatialExperiment` and `SpatialFeatureExperiment` objects
and lets users flexibly define their own `ppm`/`mppm` models with R's formula interface.

## Installation

You can install the development version of `multipointR` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mjemons/multipointR")
```

## Disclaimer

This package is still under active development, the content is therefore
subject to change. 

Parts of this code were optimised and/or generated with claude.ai. 
Whenever claude.ai made significant contributions to a code chunk 
I tried my best to acknowledge this as comments in the code.

## Contact

In case you have suggestions to `multipointR` please consider opening an
issue to this repository.