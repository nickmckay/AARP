# Log-posterior for the forward-model two-acid calibration (Stage 1)

Log-posterior for the forward-model two-acid calibration (Stage 1)

## Usage

``` r
log_posterior_age_calibration_fwd(params, precomputed, dT_yr = 1500)
```

## Arguments

- params:

  Named numeric vector: `log_Ae_Asp`, `Ea_Asp`, `x_Asp`,
  `log_sigma_Asp`, `log_Ae_Glu`, `Ea_Glu`, `x_Glu`, `log_sigma_Glu`,
  `z_rho`.

- precomputed:

  List returned by
  [`precompute_forward_inputs`](https://nickmckay.github.io/AARP/reference/precompute_forward_inputs.md).

- dT_yr:

  Integration timestep (yr). Default `1500`.

## Value

Scalar log-posterior.
