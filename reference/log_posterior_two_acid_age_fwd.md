# Log-posterior for forward-model two-acid age inversion (Stage 2)

Log-posterior for forward-model two-acid age inversion (Stage 2)

## Usage

``` r
log_posterior_two_acid_age_fwd(
  params,
  calib_params_fwd,
  DL_Asp_obs,
  DL_Glu_obs,
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

- DL_Asp_obs, DL_Glu_obs:

  Observed D/L ratios for the sample.

- temp_C:

  Representative bottom-water temperature for the sample.

- depth_cm:

  Placeholder depth (cm). Default `100`.

- dT_yr:

  Integration timestep (yr). Default `2000`.

- range_ka:

  Numeric length-2 vector passed to
  [`log_prior_age`](https://nickmckay.github.io/AARP/reference/log_prior_age.md).

## Value

Scalar log-posterior.
