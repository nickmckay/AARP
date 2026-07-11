# Predict age from D/L using the TDK (Asp) calibration curve

\$\$t = a \cdot \mathrm{atanh}(D/L)^e\$\$

## Usage

``` r
predict_age_asp(DL, a, e)
```

## Arguments

- DL:

  Numeric vector of aspartic acid D/L ratios.

- a:

  Scale parameter (ka).

- e:

  Exponent.

## Value

Numeric vector of predicted ages (ka BP).

## See also

[`predict_DL_asp`](https://nickmckay.github.io/AARP/reference/predict_DL_asp.md)
for the inverse transform.
