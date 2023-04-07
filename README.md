
# tcrpheno

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install development version of tcrpheno from Github with:

``` r
remotes::install_github("kalaga27/tcrpheno")
```

## Main Function

The main function is:

``` r
score_tcrs(data, chain)
```
arguments:

`data`: input TCR data. the tcrpheno package contains example data `tcrpheno_input` to demonstrate column names. first column is cell identifier, order of columns does not matter otherwise.

`chain`: is your TCR data paired or single-chain? `"ab"` for paired ab TCRs, `"a"` for alpha chain only, `"b"` for beta chain only

## Output

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

