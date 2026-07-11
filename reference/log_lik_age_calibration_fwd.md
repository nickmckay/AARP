# Log-likelihood for the forward-model two-acid calibration

Forward-simulates D/L for both acids at every calibration sample using a
single shared kinetic model form (Arrhenius + Bada power-law, via
[`racemize_one_depth`](https://nickmckay.github.io/AARP/reference/racemize_one_depth.md)),
each acid with its own `Ae`, `Ea`, `x`, and evaluates a correlated
bivariate-normal likelihood on the D/L residuals (paralleling
[`log_lik_age_calibration`](https://nickmckay.github.io/AARP/reference/log_lik_age_calibration.md),
but in D/L space rather than log-age space, since prediction now happens
forward from age to D/L rather than via an invertible empirical curve).

## Usage

``` r
log_lik_age_calibration_fwd(
  log_Ae_Asp,
  Ea_Asp,
  x_Asp,
  log_sigma_Asp,
  log_Ae_Glu,
  Ea_Glu,
  x_Glu,
  log_sigma_Glu,
  z_rho,
  precomputed,
  dT_yr = 1500
)
```

## Arguments

- log_Ae_Asp, Ea_Asp, x_Asp, log_sigma_Asp:

  Aspartic acid kinetic parameters and log residual SD (D/L scale).

- log_Ae_Glu, Ea_Glu, x_Glu, log_sigma_Glu:

  Glutamic acid kinetic parameters and log residual SD.

- z_rho:

  Fisher-z transform of the residual correlation.

- precomputed:

  List returned by
  [`precompute_forward_inputs`](https://nickmckay.github.io/AARP/reference/precompute_forward_inputs.md).

- dT_yr:

  Integration timestep (yr). Default `1500`.

## Value

Scalar log-likelihood.
