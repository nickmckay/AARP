# Log-prior for the forward-model two-acid calibration

`log_Ae_*` and `x_*` are weakly informative, centred loosely on the
package's existing (isoleucine) defaults
([`default_AAR_params`](https://nickmckay.github.io/AARP/reference/default_AAR_params.md)).
`Ea_Asp`/`Ea_Glu` use an **informed** prior centred in the ~25-35
kcal/mol range typical of published racemization activation energies for
calcite-hosted amino acids – deliberately more informative than the
other parameters, because the calibration dataset's narrow natural
temperature range (~3.5 deg C) cannot identify Ea on its own (see
[`vignette("two-acid-age-forward-model")`](https://nickmckay.github.io/AARP/articles/two-acid-age-forward-model.md)).

## Usage

``` r
log_prior_age_calibration_fwd(
  log_Ae_Asp,
  Ea_Asp,
  x_Asp,
  log_sigma_Asp,
  log_Ae_Glu,
  Ea_Glu,
  x_Glu,
  log_sigma_Glu,
  z_rho
)
```

## Arguments

- log_Ae_Asp, Ea_Asp, x_Asp, log_sigma_Asp:

  Aspartic acid kinetic parameters and log residual SD (D/L scale).

- log_Ae_Glu, Ea_Glu, x_Glu, log_sigma_Glu:

  Glutamic acid kinetic parameters and log residual SD.

- z_rho:

  Fisher-z transform of the residual correlation.

## Value

Scalar log-prior.
