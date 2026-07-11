# Log-likelihood for the joint two-acid age calibration

Models \\\log(t)\\ for each acid as a linear function of
\\\log(f(D/L))\\ (TDK for Asp, SPK for Glu), with residuals for the two
acids drawn from a correlated bivariate normal at each calibration
sample (both acids are measured on the same dated sample, so their
departures from the mean trend are expected to covary).

## Usage

``` r
log_lik_age_calibration(
  log_a_Asp,
  e_Asp,
  log_sigma_Asp,
  log_a_Glu,
  e_Glu,
  log_sigma_Glu,
  z_rho,
  calibration_data
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

- calibration_data:

  Data frame with columns `age_ka`, `DL_Asp`, `DL_Glu`.

## Value

Scalar log-likelihood.
