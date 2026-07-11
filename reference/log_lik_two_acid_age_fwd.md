# Log-likelihood for forward-model two-acid age inversion

Forward-simulates D/L for both acids at a candidate age – via a trivial
single-point age-depth model and a simulated placeholder temperature
history (see
[`simulate_temp_history`](https://nickmckay.github.io/AARP/reference/simulate_temp_history.md))
– and evaluates a correlated bivariate-normal likelihood on the D/L
residuals, using the shared kinetic model fit in Stage 1
([`log_posterior_age_calibration_fwd`](https://nickmckay.github.io/AARP/reference/log_posterior_age_calibration_fwd.md)).

## Usage

``` r
log_lik_two_acid_age_fwd(
  log_age,
  calib_params_fwd,
  DL_Asp_obs,
  DL_Glu_obs,
  temp_C,
  depth_cm = 100,
  dT_yr = 2000
)
```

## Arguments

- log_age:

  Candidate value of \\\log(\text{age in ka})\\.

- calib_params_fwd:

  List with elements `Ae_Asp`, `Ea_Asp`, `x_Asp`, `sigma_Asp`, `Ae_Glu`,
  `Ea_Glu`, `x_Glu`, `sigma_Glu`, `rho` (see
  [`summarize_calibration_fwd`](https://nickmckay.github.io/AARP/reference/summarize_calibration_fwd.md)).

- DL_Asp_obs, DL_Glu_obs:

  Observed D/L ratios for the sample.

- temp_C:

  Representative bottom-water temperature for the sample.

- depth_cm:

  Placeholder depth (cm) used only to define the trivial age-depth
  model; the result does not depend on its value. Default `100`.

- dT_yr:

  Integration timestep (yr). Default `2000`.

## Value

Scalar log-likelihood.
