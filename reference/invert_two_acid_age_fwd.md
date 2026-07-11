# Run the forward-model two-acid age inversion for one sample

Convenience wrapper around
[`run_mcmc`](https://nickmckay.github.io/AARP/reference/run_mcmc.md) for
[`log_posterior_two_acid_age_fwd`](https://nickmckay.github.io/AARP/reference/log_posterior_two_acid_age_fwd.md).
Produces the same single-parameter (`log_age`) chain shape as
[`invert_two_acid_age`](https://nickmckay.github.io/AARP/reference/invert_two_acid_age.md),
so
[`summarize_age_posterior`](https://nickmckay.github.io/AARP/reference/summarize_age_posterior.md)
and
[`plot_age_posterior`](https://nickmckay.github.io/AARP/reference/plot_age_posterior.md)
apply unchanged.

## Usage

``` r
invert_two_acid_age_fwd(
  DL_Asp_obs,
  DL_Glu_obs,
  calib_params_fwd,
  temp_C,
  depth_cm = 100,
  n_iter = 20000,
  proposal_sd = 0.3,
  dT_yr = 2000,
  range_ka = c(0.5, 3000)
)
```

## Arguments

- DL_Asp_obs, DL_Glu_obs:

  Observed D/L ratios for the sample.

- calib_params_fwd:

  List of calibration constants, see
  [`log_lik_two_acid_age_fwd`](https://nickmckay.github.io/AARP/reference/log_lik_two_acid_age_fwd.md).

- temp_C:

  Representative bottom-water temperature for the sample.

- depth_cm:

  Placeholder depth (cm). Default `100`.

- n_iter:

  Integer. Total MCMC iterations. Default `20000`.

- proposal_sd:

  Proposal SD for `log_age`. Default `0.3`.

- dT_yr:

  Integration timestep (yr). Default `2000`.

- range_ka:

  Numeric length-2 vector passed to
  [`log_prior_age`](https://nickmckay.github.io/AARP/reference/log_prior_age.md).

## Value

List returned by
[`run_mcmc`](https://nickmckay.github.io/AARP/reference/run_mcmc.md).
