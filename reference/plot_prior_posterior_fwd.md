# Plot prior vs. posterior densities for all Stage 1 forward-model parameters

Overlays each parameter's prior density (as specified in
[`log_prior_age_calibration_fwd`](https://nickmckay.github.io/AARP/reference/log_prior_age_calibration_fwd.md))
with a posterior density estimated from the post-burn-in MCMC samples. A
posterior much narrower than its prior means the data are informative
for that parameter; a posterior that closely tracks the prior (as for
`Ea_Asp`/`Ea_Glu`, see
[`plot_ae_ea_posterior`](https://nickmckay.github.io/AARP/reference/plot_ae_ea_posterior.md))
means the data alone cannot move it far from the assumed prior belief.

## Usage

``` r
plot_prior_posterior_fwd(mcmc_out, burnin = 1000)
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
