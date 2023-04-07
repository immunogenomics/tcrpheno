
# tcrpheno

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install development version of tcrpheno from Github with:

``` r
remotes::install_github("kalaga27/tcrpheno")
```

## Example

This is a basic example:

``` r
library(tcrpheno)
head(tcrpheno_input)
result = score_tcrs(tcrpheno_input, chain="ab")
head(result)
```

