<!-- badges: start -->
[![R-CMD-check](https://github.com/andrewherren/bayeslm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/andrewherren/bayeslm/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# bayeslm

## Description

Efficient sampling for Gaussian linear regression with arbitrary priors. This package implements Bayesian linear regression using elliptical slice sampler, which allows easily usage of arbitrary priors. This package is parallelized by RcppParallel.

## Installation

```
install.packages("devtools")
library(devtools)
install_github("JingyuHe/bayeslm")
```

## Reference

The method underlying this package is described in "Efficient sampling for Gaussian linear regression with arbitrary priors" (Hahn, He, and Lopes 2019) which was [published in the Journal of Computational and Graphical Statistics](https://www.tandfonline.com/doi/abs/10.1080/10618600.2018.1482762).

An open-access version of the paper is [available on Arxiv](https://arxiv.org/abs/1806.05738).
