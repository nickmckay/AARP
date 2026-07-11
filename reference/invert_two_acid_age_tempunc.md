# Run the two-acid age inversion with temperature uncertainty folded in

Run the two-acid age inversion with temperature uncertainty folded in

## Usage

``` r
invert_two_acid_age_tempunc(
  DL_Asp_obs,
  DL_Glu_obs,
  calib_params_fwd,
  temp_C_center,
  depth_cm = 100,
  n_iter = 20000,
  proposal_sd = c(log_age = 0.3, temp_offset = 0.3),
  dT_yr = 2000,
  range_ka = c(0.5, 3000),
  temp_sd = 1
)
```

## Arguments

- DL_Asp_obs, DL_Glu_obs:

  Observed D/L ratios for the sample.

- calib_params_fwd:

  List of calibration constants, see
  [`log_lik_two_acid_age_fwd`](https://nickmckay.github.io/AARP/reference/log_lik_two_acid_age_fwd.md).

- temp_C_center:

  Central/recorded bottom-water temperature (deg C).

- depth_cm:

  Placeholder depth (cm). Default `100`.

- n_iter:

  Integer. Total MCMC iterations. Default `20000`.

- proposal_sd:

  Named numeric vector of proposal SDs for `log_age` and `temp_offset`.
  Default `c(log_age = 0.3, temp_offset = 0.3)`.

- dT_yr:

  Integration timestep (yr). Default `2000`.

- range_ka:

  Numeric length-2 vector passed to
  [`log_prior_age`](https://nickmckay.github.io/AARP/reference/log_prior_age.md).

- temp_sd:

  Prior SD on `temp_offset` (deg C). Default `1`.

## Value

List returned by
[`run_mcmc`](https://nickmckay.github.io/AARP/reference/run_mcmc.md).
