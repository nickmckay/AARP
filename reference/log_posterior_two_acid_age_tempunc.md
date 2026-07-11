# Log-posterior for two-acid age inversion with temperature uncertainty

Log-posterior for two-acid age inversion with temperature uncertainty

## Usage

``` r
log_posterior_two_acid_age_tempunc(
  params,
  calib_params_fwd,
  DL_Asp_obs,
  DL_Glu_obs,
  temp_C_center,
  depth_cm = 100,
  dT_yr = 2000,
  range_ka = c(0.5, 3000),
  temp_sd = 1
)
```

## Arguments

- params:

  Named numeric vector with elements `log_age`, `temp_offset`.

- calib_params_fwd:

  List of calibration constants, see
  [`log_lik_two_acid_age_fwd`](https://nickmckay.github.io/AARP/reference/log_lik_two_acid_age_fwd.md).

- DL_Asp_obs, DL_Glu_obs:

  Observed D/L ratios for the sample.

- temp_C_center:

  Central/recorded bottom-water temperature (deg C).

- depth_cm, dT_yr:

  As in
  [`log_lik_two_acid_age_fwd`](https://nickmckay.github.io/AARP/reference/log_lik_two_acid_age_fwd.md).

- range_ka:

  Numeric length-2 vector passed to
  [`log_prior_age`](https://nickmckay.github.io/AARP/reference/log_prior_age.md).

- temp_sd:

  Prior SD on `temp_offset` (deg C). Default `1`, reflecting a plausible
  cold/narrow-but-not-precisely-known range.

## Value

Scalar log-posterior.
