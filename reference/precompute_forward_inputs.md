# Precompute per-sample forward-model inputs for calibration

Builds the (depth, age-model, temperature-history) triple for each
calibration sample once, so the expensive parts of the forward model
(age-depth interpolation, temperature-history construction) are not
repeated on every MCMC iteration – only the kinetic parameters change
between iterations.

## Usage

``` r
precompute_forward_inputs(
  calibration_data,
  amplitude = 1,
  period_ka = 100,
  resolution_ka = 5
)
```

## Arguments

- calibration_data:

  Data frame with columns `study`, `core`, `depth_mbsf`, `age_ka`,
  `temp_C`, `DL_Asp`, `DL_Glu`. Rows with missing `temp_C` are dropped.

- amplitude, period_ka, resolution_ka:

  Passed to
  [`simulate_temp_history`](https://nickmckay.github.io/AARP/reference/simulate_temp_history.md).

## Value

A list with elements `depth_cm`, `age_model` (list), `temp_model`
(list), `DL_Asp`, `DL_Glu`, one entry/row per retained sample.
