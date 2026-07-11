# Log-posterior for two-acid age inversion (Stage 2)

Log-posterior for two-acid age inversion (Stage 2)

## Usage

``` r
log_posterior_two_acid_age(
  params,
  calib_params,
  DL_Asp_obs,
  DL_Glu_obs,
  range_ka = c(0.5, 3000)
)
```

## Arguments

- params:

  Named numeric vector with element `log_age`.

- calib_params:

  List of calibration constants, see
  [`log_lik_two_acid_age`](https://nickmckay.github.io/AARP/reference/log_lik_two_acid_age.md).

- DL_Asp_obs, DL_Glu_obs:

  Observed D/L ratios for the sample.

- range_ka:

  Numeric length-2 vector passed to
  [`log_prior_age`](https://nickmckay.github.io/AARP/reference/log_prior_age.md).

## Value

Scalar log-posterior.
