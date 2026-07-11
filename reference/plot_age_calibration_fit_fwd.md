# Plot forward-model predicted vs. observed D/L for the calibration set

Plot forward-model predicted vs. observed D/L for the calibration set

## Usage

``` r
plot_age_calibration_fit_fwd(
  mcmc_out,
  precomputed,
  burnin = 1000,
  dT_yr = 1500
)
```

## Arguments

- mcmc_out:

  List returned by
  [`run_mcmc`](https://nickmckay.github.io/AARP/reference/run_mcmc.md)
  using
  [`log_posterior_age_calibration_fwd`](https://nickmckay.github.io/AARP/reference/log_posterior_age_calibration_fwd.md).

- precomputed:

  List returned by
  [`precompute_forward_inputs`](https://nickmckay.github.io/AARP/reference/precompute_forward_inputs.md).

- burnin:

  Integer. Burn-in iterations to discard. Default `1000`.

- dT_yr:

  Integration timestep (yr). Default `1500`.

## Value

A ggplot object.
