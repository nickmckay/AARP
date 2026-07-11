# Bivariate normal log-density with zero means

Bivariate normal log-density with zero means

## Usage

``` r
.dbvnorm_log(x1, x2, sigma1, sigma2, rho)
```

## Arguments

- x1, x2:

  Numeric vectors of residuals (observed minus predicted).

- sigma1, sigma2:

  Marginal standard deviations.

- rho:

  Correlation coefficient, in \\(-1, 1)\\.

## Value

Numeric vector of log-densities (summed if used with
[`sum()`](https://rdrr.io/r/base/sum.html)).
