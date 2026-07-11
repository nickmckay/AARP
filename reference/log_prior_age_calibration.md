# Log-prior for the joint two-acid age calibration

Weakly informative priors centred on literature-scale values (Kaufman et
al. 2013-type Np calibrations), wide enough that a calibration dataset
of ~100+ samples dominates the posterior.

## Usage

``` r
log_prior_age_calibration(
  log_a_Asp,
  e_Asp,
  log_sigma_Asp,
  log_a_Glu,
  e_Glu,
  log_sigma_Glu,
  z_rho
)
```

## Arguments

- log_a_Asp, e_Asp, log_sigma_Asp:

  TDK (Asp) calibration parameters: log scale factor, exponent, log
  residual SD (log-age scale).

- log_a_Glu, e_Glu, log_sigma_Glu:

  SPK (Glu) calibration parameters.

- z_rho:

  Fisher-z transform of the residual correlation (\\\rho =
  \tanh(z\_\rho)\\).

## Value

Scalar log-prior.
