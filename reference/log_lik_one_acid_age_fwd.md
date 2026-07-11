# Log-likelihood for single-acid forward-model age inversion

Same forward model as
[`log_lik_two_acid_age_fwd`](https://nickmckay.github.io/AARP/reference/log_lik_two_acid_age_fwd.md),
but using only one acid's observation and its own (univariate) residual
SD – useful for quantifying how much a second acid narrows the age
estimate.

## Usage

``` r
log_lik_one_acid_age_fwd(
  log_age,
  calib_params_fwd,
  DL_obs,
  acid = c("Asp", "Glu"),
  temp_C,
  depth_cm = 100,
  dT_yr = 2000
)
```

## Arguments

- log_age:

  Candidate value of \\\log(\text{age in ka})\\.

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

## Value

Scalar log-likelihood.
