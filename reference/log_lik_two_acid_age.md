# Log-likelihood for two-acid age inversion

Converts a sample's observed Asp and Glu D/L into two correlated point
estimates of \\\log(\text{age})\\ via the Stage 1 calibration curves
([`predict_age_asp`](https://nickmckay.github.io/AARP/reference/predict_age_asp.md),
[`predict_age_glu`](https://nickmckay.github.io/AARP/reference/predict_age_glu.md)),
then evaluates how consistent a candidate age is with both, using the
same residual covariance (`sigma_Asp`, `sigma_Glu`, `rho`) estimated
during calibration. This is the Bayesian analogue of Kaufman's GLS
combination of Asp and Glu ages: instead of algebraically combining two
point estimates and variances (with an ad hoc floor to keep GLS weights
positive), the two acids simply contribute correlated evidence to one
joint likelihood, which combines coherently for any value of `rho`.

## Usage

``` r
log_lik_two_acid_age(log_age, calib_params, DL_Asp_obs, DL_Glu_obs)
```

## Arguments

- log_age:

  Candidate value of \\\log(\text{age in ka})\\.

- calib_params:

  List with elements `a_Asp`, `e_Asp`, `sigma_Asp`, `a_Glu`, `e_Glu`,
  `sigma_Glu`, `rho` (see
  [`summarize_calibration`](https://nickmckay.github.io/AARP/reference/summarize_calibration.md)).

- DL_Asp_obs, DL_Glu_obs:

  Observed D/L ratios for the sample.

## Value

Scalar log-likelihood.
