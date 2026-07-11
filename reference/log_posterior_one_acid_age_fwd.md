# Log-posterior for single-acid forward-model age inversion

Log-posterior for single-acid forward-model age inversion

## Usage

``` r
log_posterior_one_acid_age_fwd(
  params,
  calib_params_fwd,
  DL_obs,
  acid = c("Asp", "Glu"),
  temp_C,
  depth_cm = 100,
  dT_yr = 2000,
  range_ka = c(0.5, 3000)
)
```

## Arguments

- params:

  Named numeric vector with element `log_age`.

- calib_params_fwd:

  List of calibration constants, see
  [`log_lik_two_acid_age_fwd`](https://nickmckay.github.io/AARP/reference/log_lik_two_acid_age_fwd.md).

- DL_obs:

  Observed D/L ratio for the chosen acid.

- acid:

  Character, `"Asp"` or `"Glu"`.

- temp_C, depth_cm, dT_yr:

  As in
  [`log_lik_two_acid_age_fwd`](https://nickmckay.github.io/AARP/reference/log_lik_two_acid_age_fwd.md).

- range_ka:

  Numeric length-2 vector passed to
  [`log_prior_age`](https://nickmckay.github.io/AARP/reference/log_prior_age.md).

## Value

Scalar log-posterior.
