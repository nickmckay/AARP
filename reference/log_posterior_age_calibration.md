# Log-posterior for the joint two-acid age calibration (Stage 1)

Log-posterior for the joint two-acid age calibration (Stage 1)

## Usage

``` r
log_posterior_age_calibration(params, calibration_data)
```

## Arguments

- params:

  Named numeric vector: `log_a_Asp`, `e_Asp`, `log_sigma_Asp`,
  `log_a_Glu`, `e_Glu`, `log_sigma_Glu`, `z_rho`.

- calibration_data:

  Data frame with columns `age_ka`, `DL_Asp`, `DL_Glu`.

## Value

Scalar log-posterior (possibly `-Inf`).
