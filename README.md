
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ATbounds

<!-- badges: start -->

[![](https://cranlogs.r-pkg.org/badges/ATbounds)](https://cran.r-project.org/package=ATbounds)
[![codecov](https://codecov.io/gh/ATbounds/ATbounds-r/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ATbounds/ATbounds-r)
<!-- badges: end -->

***ATbounds*** is an R package that provides estimation and inference
methods for bounding average treatment effects (on the treated) that are
valid under an unconfoundedness assumption. The bounds are designed to
be robust in challenging situations, for example, when the conditioning
variables take on a large number of different values in the observed
sample, or when the overlap condition is violated. This robustness is
achieved by only using limited “pooling” of information across
observations.

For more details, see Lee and Weidner, “Bounding Treatment Effects by
Pooling Limited Information across Observations,” forthcoming at the
*Journal of Econometrics* and available at [arXiv:2111.05243
\[econ.EM\]](https://arxiv.org/abs/2111.05243).

## Installation

You can install the released version of ***ATbounds*** from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("ATbounds")
```

Alternatively, you can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools") if devtools is not installed
devtools::install_github("https://github.com/ATbounds/ATbounds-r")
```

## Documentation and Examples

The package vignette provides worked examples using the
NSW/Dehejia-Wahba, right heart catheterization, and electronic fetal
monitoring datasets. After installing the package, the vignette can be
opened in R with:

``` r
vignette("ATbounds_vignette", package = "ATbounds")
```

## Replication Materials

The companion repository
[ATbounds/ATbounds-replication](https://github.com/ATbounds/ATbounds-replication)
contains the replication materials for Lee and Weidner, “Bounding
Treatment Effects by Pooling Limited Information across Observations,”
forthcoming at the *Journal of Econometrics*.

The replication repository contains the scripts and generated outputs
used to reproduce the paper’s numerical results. It is separate from
this package repository: the vignette illustrates how to use the
***ATbounds*** package, while the replication repository reproduces the
paper tables and main-text numerical results.

## Reference

Lee, S. and Weidner, M. (forthcoming). “Bounding Treatment Effects by
Pooling Limited Information across Observations.” *Journal of
Econometrics*. [arXiv:2111.05243
\[econ.EM\]](https://arxiv.org/abs/2111.05243).
