
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MultiNeSS

<!-- badges: start -->

<!-- badges: end -->

implements model fitting and simulation for Gaussian and logistic inner
product MultiNeSS models for multiplex networks. The package uses a
convex fitting algorithm with fully adaptive parameter tuning, including
options for edge cross-validation. For more details see [MacDonald et
al., (2020)](https://arxiv.org/abs/2012.14409).

## Installation

You can install the development version of  from GitHub using

``` r
devtools::install_github("peterwmacd/multiness")
```

## Example

includes an example agricultural trade dataset which is studied in
[MacDonald et al., (2020)](https://arxiv.org/abs/2012.14409). It is easy
to import and to fit a Gaussian MultiNeSS model with adaptive tuning.

``` r
library(multiness)

# import data
data(agri_trade)
# log transformation for edge weights
A <- log(1+agri_trade)

# model fit
fit <- multiness_fit(A=A,model="gaussian",self_loops=FALSE,
                     tuning="adaptive",tuning_opts=list(penalty_const=3),
                     optim_opts=list(max_rank=100))

# inspect fitted latent space dimensions
# common latent space
fit$d1
#> [1] 23
# individual latent spaces
fit$d2
#>  [1] 3 3 4 1 2 4 3 4 8 7 3 6 3
```
