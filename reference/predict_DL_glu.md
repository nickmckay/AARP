# Predict D/L from age using the SPK (Glu) calibration curve

Inverse of
[`predict_age_glu`](https://nickmckay.github.io/AARP/reference/predict_age_glu.md):
\\D/L = (t/a)^{1/e}\\.

## Usage

``` r
predict_DL_glu(age_ka, a, e)
```

## Arguments

- age_ka:

  Numeric vector of ages (ka BP).

- a:

  Scale parameter (ka).

- e:

  Exponent.

## Value

Numeric vector of predicted D/L (Glu) ratios.
