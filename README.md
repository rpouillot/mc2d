[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/mc2d)](https://cran.r-project.org/package=mc2d)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/mc2d)](https://cran.r-project.org/package=mc2d)

# Tools for Two-Dimensional Monte-Carlo Simulations

## The package 

`mc2d` provides a complete framework to build and study Two-Dimensional Monte-Carlo simulations, aka Second-Order Monte-Carlo simulations. It also includes various distributions frequently used in the risk assessment domain (pert, triangular, Bernoulli, empirical discrete and continuous, beta subjective, 
Minimum Quantile Information Distribution, ...).



## Getting it

The stable version of `mc2d` can be installed from CRAN using:
```r
install.packages("mc2d")
library("mc2d")
```

The development version of `mc2d` can be installed from GitHub (`devtools` needed):

```r
if (!requireNamespace("devtools", quietly = TRUE))
   install.packages("devtools")
devtools::install_github("ericcheny/mc2d-update")
library("mc2d")
```

Check the NEWS [here](https://github.com/rpouillot/mc2d/blob/main/inst/NEWS).

## Documentation

See the manual and the vignettes distributed with the package.

## Issues

Issues can be reported on https://github.com/rpouillot/mc2d/issues.

or directly to the maintainer Régis Pouillot: rpouillot@yahoo.fr


## Citations

If you use `mc2d`, please cite:

R. Pouillot, M.-L. Delignette-Muller (2010), [Evaluating variability and uncertainty in microbial
  quantitative risk assessment using two R packages](www.doi.org/10.1016/j.ijfoodmicro.2010.07.011). International Journal of Food Microbiology.
  142(3):330-40  
