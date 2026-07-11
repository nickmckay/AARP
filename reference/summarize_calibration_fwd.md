# Extract posterior-mean forward-model calibration parameters

Extract posterior-mean forward-model calibration parameters

## Usage

``` r
summarize_calibration_fwd(mcmc_out, burnin = 1000)
```

## Arguments

- mcmc_out:

  List returned by
  [`run_mcmc`](https://nickmckay.github.io/AARP/reference/run_mcmc.md)
  using
  [`log_posterior_age_calibration_fwd`](https://nickmckay.github.io/AARP/reference/log_posterior_age_calibration_fwd.md).

- burnin:

  Integer. Burn-in iterations to discard. Default `1000`.

## Value

A named list with elements `Ae_Asp`, `Ea_Asp`, `x_Asp`, `sigma_Asp`,
`Ae_Glu`, `Ea_Glu`, `x_Glu`, `sigma_Glu`, `rho`.
