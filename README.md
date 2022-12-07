
# fleissci <img src="man/figures/logo.png" align="right" width="411" height="321"/>

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/fleissci.png)](https://cran.r-project.org/package=fleissci)
[![R-CMD-check](https://github.com/JonasMoss/fleissci/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JonasMoss/fleissci/actions/workflows/R-CMD-check.yaml)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Codecov test
coverage](https://codecov.io/gh/JonasMoss/fleissci/branch/main/graph/badge.svg)](https://app.codecov.io/gh/JonasMoss/fleissci?branch=main)

An `R` package for calculating and doing inference with the
quadratically weighted Fleiss’ kappa, Cohen’s kappa, and Conger’s kappa.

**Note:** I justed started development on this package!

## Installation

The package is not available on `CRAN` yet, so use the following command
from inside `R`:

``` r
# install.packages("remotes")
remotes::install_github("JonasMoss/fleissci")
```

## Usage

Call the `library` function and load the data of Zapf et al. (2016):

``` r
library("fleissci")
head(dat.zapf2016)
#>   Rater A Rater B Rater C Rater D
#> 1       5       5       4       5
#> 2       1       1       1       1
#> 3       5       5       5       5
#> 4       1       3       3       3
#> 5       5       5       5       5
#> 6       1       1       1       1
```

Then calculate an asymptotically distribution-free confidence interval
for ![\kappa](https://latex.codecogs.com/svg.latex?%5Ckappa "\kappa"),

``` r
fleissci(dat.zapf2016)
#> Call: fleissci(x = dat.zapf2016)
#> 
#> 95% confidence interval (n = 50).
#>     0.025     0.975 
#> 0.8429321 0.9538451 
#> 
#> Sample estimates.
#>     kappa        sd 
#> 0.8983886 0.1931729
```

You can also calculate confidence intervals for Conger’s kappa (Cohen’s
kappa.)

``` r
congerci(dat.zapf2016)
#> Call: congerci(x = dat.zapf2016)
#> 
#> 95% confidence interval (n = 50).
#>     0.025     0.975 
#> 0.8430854 0.9538547 
#> 
#> Sample estimates.
#>     kappa        sd 
#> 0.8984700 0.1929226
```

## Supported techniques

`fleissci` supports three basic asymptotic confidence interval
constructios. The asymptotically distribution-free interval, the
pseudo-elliptical interval, and the normal method.

| Method       | Description                                                                                                                                                                                                                                                                                        |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `adf`        | The asymptotic distribution free method. The method is asymptotically correct, but has poor small-sample performance.                                                                                                                                                                              |
| `elliptical` | The elliptical or pseudo-elliptical kurtosis correction. Uses the unbiased sample estimator of the common kurtosis (Joanes, 1998). Has better small-sample performance than `adf` and `normal` if the kurtosis is large and ![n](https://latex.codecogs.com/svg.latex?n "n") is small.             |
| `normal`     | Assumes normality of ![X](https://latex.codecogs.com/svg.latex?X "X"). This method is not recommended since it yields too short confidence intervals when the excess kurtosis of ![X](https://latex.codecogs.com/svg.latex?X "X") is larger than ![0](https://latex.codecogs.com/svg.latex?0 "0"). |

In addition, you may transform the intervals using one of four
transforms:

1.  The [Fisher
    transform](https://en.wikipedia.org/wiki/Fisher_transformation), or
    ![\kappa\mapsto \operatorname{artanh}(\kappa)](https://latex.codecogs.com/svg.latex?%5Ckappa%5Cmapsto%20%5Coperatorname%7Bartanh%7D%28%5Ckappa%29 "\kappa\mapsto \operatorname{artanh}(\kappa)").
    Famously used in inference for the correlation coefficient.
2.  The ![\log](https://latex.codecogs.com/svg.latex?%5Clog "\log")
    transform, where
    ![\kappa \mapsto \log(1-\kappa)](https://latex.codecogs.com/svg.latex?%5Ckappa%20%5Cmapsto%20%5Clog%281-%5Ckappa%29 "\kappa \mapsto \log(1-\kappa)").
    This is an asymptotic pivot under the elliptical model with parallel
    items.
3.  The identity transform. The default option.
4.  The
    [![\arcsin](https://latex.codecogs.com/svg.latex?%5Carcsin "\arcsin")
    transform](https://en.wikipedia.org/wiki/Inverse_trigonometric_functions).
    This transform might fail when
    ![n](https://latex.codecogs.com/svg.latex?n "n") is small, as
    negative values for
    ![\hat{\kappa}](https://latex.codecogs.com/svg.latex?%5Chat%7B%5Ckappa%7D "\hat{\kappa}")
    is possible, but
    ![\arcsin](https://latex.codecogs.com/svg.latex?%5Carcsin "\arcsin")
    do not accept them,

The option `bootstrap` does studentized bootstrapping Efron, B. (1987)
with `n_reps` repetitions. If `bootstrap = FALSE`, an ordinary normal
approximation will be used. The studentized bootstrap intervals are is a
second-order correct, so its confidence intervals will be better than
the normal approximation when
![n](https://latex.codecogs.com/svg.latex?n "n") is sufficiently large.

## Similar software

There are several `R` packages that calculate agreement coefficients,
most notably `irrCAC`.

## How to Contribute or Get Help

If you encounter a bug, have a feature request or need some help, open a
[Github issue](https://github.com/JonasMoss/fleissci/issues). Create a
pull requests to contribute.

## References

- Zapf, A., Castell, S., Morawietz, L. et al. Measuring inter-rater
  reliability for nominal data – which coefficients and confidence
  intervals are appropriate?. BMC Med Res Methodol 16, 93 (2016).
  https://doi.org/10.1186/s12874-016-0200-9
