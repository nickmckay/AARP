# Run the two-acid Bayesian age inversion for one sample

Convenience wrapper around
[`run_mcmc`](https://nickmckay.github.io/AARP/reference/run_mcmc.md) for
[`log_posterior_two_acid_age`](https://nickmckay.github.io/AARP/reference/log_posterior_two_acid_age.md).
Initialises the chain at the average of the two single-acid point-age
estimates.

## Usage

``` r
invert_two_acid_age(
  DL_Asp_obs,
  DL_Glu_obs,
  calib_params,
  n_iter = 20000,
  proposal_sd = 0.3,
  range_ka = c(0.5, 3000)
)
```

## Arguments

- DL_Asp_obs, DL_Glu_obs:

  Observed D/L ratios for the sample.

- calib_params:

  List of calibration constants, see
  [`log_lik_two_acid_age`](https://nickmckay.github.io/AARP/reference/log_lik_two_acid_age.md).

- n_iter:

  Integer. Total MCMC iterations. Default `20000`.

- proposal_sd:

  Proposal SD for `log_age`. Default `0.3`.

- range_ka:

  Numeric length-2 vector passed to
  [`log_prior_age`](https://nickmckay.github.io/AARP/reference/log_prior_age.md).

## Value

List returned by
[`run_mcmc`](https://nickmckay.github.io/AARP/reference/run_mcmc.md).
