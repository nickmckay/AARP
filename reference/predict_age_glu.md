# Predict age from D/L using the SPK (Glu) calibration curve

\$\$t = a \cdot (D/L)^e\$\$

## Usage

``` r
predict_age_glu(DL, a, e)
```

## Arguments

- DL:

  Numeric vector of glutamic acid D/L ratios.

- a:

  Scale parameter (ka).

- e:

  Exponent.

## Value

Numeric vector of predicted ages (ka BP).

## See also

[`predict_DL_glu`](https://nickmckay.github.io/AARP/reference/predict_DL_glu.md)
for the inverse transform.
