# Plot the Ae-Ea posterior compensation ridge

Amino acid racemization rate constants show a classic Arrhenius
"compensation effect": with a narrow natural temperature range, `Ae` and
`Ea` become strongly (positively) correlated in the posterior, because
many `(Ae, Ea)` combinations give nearly the same rate constant over the
small span of temperatures actually observed. This plot makes that
identifiability limitation visible directly.

## Usage

``` r
plot_ae_ea_posterior(mcmc_out, burnin = 1000)
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

A ggplot object.
