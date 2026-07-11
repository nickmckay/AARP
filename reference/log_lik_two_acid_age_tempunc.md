# Log-likelihood for two-acid age inversion with an uncertain temperature

As
[`log_lik_two_acid_age_fwd`](https://nickmckay.github.io/AARP/reference/log_lik_two_acid_age_fwd.md),
but the temperature history is anchored on `temp_C_center + temp_offset`
rather than a fixed value – `temp_offset` is a free parameter
representing genuine uncertainty in the sample's effective bottom-water
temperature (we know it's cold and in a narrow range, but not precisely
where in that range). Sampling `temp_offset` jointly with `log_age`
marginalizes over plausible temperature histories instead of
conditioning on one guess.

## Usage

``` r
log_lik_two_acid_age_tempunc(
  log_age,
  temp_offset,
  calib_params_fwd,
  DL_Asp_obs,
  DL_Glu_obs,
  temp_C_center,
  depth_cm = 100,
  dT_yr = 2000
)
```

## Arguments

- log_age:

  Candidate value of \\\log(\text{age in ka})\\.

- temp_offset:

  Candidate temperature offset from `temp_C_center` (deg C).

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

## Value

Scalar log-likelihood.
