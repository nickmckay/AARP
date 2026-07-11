# Predict D/L from age using the TDK (Asp) calibration curve

Inverse of
[`predict_age_asp`](https://nickmckay.github.io/AARP/reference/predict_age_asp.md):
\\D/L = \tanh((t/a)^{1/e})\\.

## Usage

``` r
predict_DL_asp(age_ka, a, e)
```

## Arguments

- age_ka:

  Numeric vector of ages (ka BP).

- a:

  Scale parameter (ka).

- e:

  Exponent.

## Value

Numeric vector of predicted D/L (Asp) ratios.
