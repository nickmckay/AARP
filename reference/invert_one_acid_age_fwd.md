# Run the single-acid forward-model age inversion for one sample

Run the single-acid forward-model age inversion for one sample

## Usage

``` r
invert_one_acid_age_fwd(
  DL_obs,
  acid = c("Asp", "Glu"),
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

- DL_obs:

  Observed D/L ratio for the chosen acid.

- acid:

  Character, `"Asp"` or `"Glu"`.

- calib_params_fwd:

  List of calibration constants, see
  [`log_lik_two_acid_age_fwd`](https://nickmckay.github.io/AARP/reference/log_lik_two_acid_age_fwd.md).

- temp_C, depth_cm, dT_yr:

  As in
  [`log_lik_two_acid_age_fwd`](https://nickmckay.github.io/AARP/reference/log_lik_two_acid_age_fwd.md).

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
