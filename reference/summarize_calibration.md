# Extract posterior-mean calibration parameters from a Stage 1 MCMC run

Convenience function that summarises a
[`run_mcmc`](https://nickmckay.github.io/AARP/reference/run_mcmc.md)
calibration chain into a single named list of calibration constants
(posterior means after burn-in), on the natural parameter scale, ready
to pass to
[`log_posterior_two_acid_age`](https://nickmckay.github.io/AARP/reference/log_posterior_two_acid_age.md)
for Stage 2 inversion.

## Usage

``` r
summarize_calibration(mcmc_out, burnin = 1000)
```

## Arguments

- mcmc_out:

  List returned by
  [`run_mcmc`](https://nickmckay.github.io/AARP/reference/run_mcmc.md)
  using
  [`log_posterior_age_calibration`](https://nickmckay.github.io/AARP/reference/log_posterior_age_calibration.md).

- burnin:

  Integer. Burn-in iterations to discard. Default `1000`.

## Value

A named list with elements `a_Asp`, `e_Asp`, `sigma_Asp`, `a_Glu`,
`e_Glu`, `sigma_Glu`, `rho`.
