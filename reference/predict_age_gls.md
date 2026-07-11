# GLS-combined two-acid age with softplus weight smoothing

Reference implementation of Darrell Kaufman's spreadsheet age
calculator: independent TDK (Asp) and SPK (Glu) age predictions are
combined by generalized least squares, with the (possibly negative) raw
GLS weights passed through a softplus floor before combination. See
[`gls_calibration_constants`](https://nickmckay.github.io/AARP/reference/gls_calibration_constants.md)
for the underlying calibration constants.

## Usage

``` r
predict_age_gls(
  DL_Asp,
  SD_Asp,
  N_Asp,
  DL_Glu,
  SD_Glu,
  N_Glu,
  constants = gls_calibration_constants
)
```

## Arguments

- DL_Asp, DL_Glu:

  Observed D/L ratios (sample means).

- SD_Asp, SD_Glu:

  Standard deviation across subsamples.

- N_Asp, N_Glu:

  Number of subsamples used for the mean.

- constants:

  List of calibration constants, see
  [`gls_calibration_constants`](https://nickmckay.github.io/AARP/reference/gls_calibration_constants.md).
  Default uses the package's built-in Np constants.

## Value

A one-row data frame with columns `t_Asp`, `t_Glu`, `V_Asp`, `V_Glu`,
`C`, `w_Asp`, `w_Glu`, `D`, `V_int`, `V_excess`, `V_total`,
`sigma_comb`, `t_comb`, `lo50`, `hi50`, `lo66`, `hi66`, `lo90`, `hi90`.

## Examples

``` r
# Reproduces "Example 1" from the spreadsheet calculator
predict_age_gls(DL_Asp = 0.350, SD_Asp = 0.01, N_Asp = 5,
                DL_Glu = 0.153, SD_Glu = 0.015, N_Glu = 5)
#>      t_Asp   t_Glu    V_Asp    V_Glu         C      w_Asp       w_Glu         D
#> 1 178.7687 153.909 0.168024 0.232258 0.1629648 0.06929328 0.007239315 0.0765326
#>       V_int    V_excess   V_total sigma_comb   t_comb     lo50     hi50
#> 1 0.1629031 0.001920084 0.1648232  0.4059842 176.2546 133.9283 231.9575
#>       lo66     hi66     lo90     hi90
#> 1 119.4717 260.0254 89.94501 345.3852
```
