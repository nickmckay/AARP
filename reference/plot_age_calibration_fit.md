# Plot fitted Asp and Glu calibration curves against calibration data

Posterior-mean TDK (Asp) and SPK (Glu) curves with 66% prediction bands
(using posterior residual SDs), overlaid on the calibration D/L-age
data.

## Usage

``` r
plot_age_calibration_fit(mcmc_out, calibration_data, burnin = 1000)
```

## Arguments

- mcmc_out:

  List returned by
  [`run_mcmc`](https://nickmckay.github.io/AARP/reference/run_mcmc.md)
  using
  [`log_posterior_age_calibration`](https://nickmckay.github.io/AARP/reference/log_posterior_age_calibration.md).

- calibration_data:

  Data frame with columns `age_ka`, `DL_Asp`, `DL_Glu`.

- burnin:

  Integer. Burn-in iterations to discard. Default `1000`.

## Value

A ggplot object.
