
# tcrpheno

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install development version of tcrpheno from Github with:

``` r
remotes::install_github("kalaga27/tcrpheno")
```

score_tcrs computes four TCR scores:

-TCRinnate: higher-scoring T cells are more likely to reach an innate-like PLZFhigh fate (MAIT and NKT TCRs score quite high)

-TCR.8: higher-scoring T cells are more likely to reach a CD8 (vs. CD4) fate

-TCRmem: higher-scoring T cells are more likely to reach a memory (vs. naive) fate

-TCRreg: higher-scoring T cells are more likely to acquire a Treg fate


## Example

This is a basic example:

``` r
library(tcrpheno)
head(tcrpheno_input)
result = score_tcrs(tcrpheno_input, chain="ab")
head(result)
```

