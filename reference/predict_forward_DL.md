# Forward-model D/L for a set of precomputed calibration samples

Runs
[`racemize_one_depth`](https://nickmckay.github.io/AARP/reference/racemize_one_depth.md)
once per sample for a given set of kinetic parameters, reusing the
precomputed age-depth models and temperature histories from
[`precompute_forward_inputs`](https://nickmckay.github.io/AARP/reference/precompute_forward_inputs.md).

## Usage

``` r
predict_forward_DL(precomputed, AAR_params, dT_yr = 1500)
```

## Arguments

- precomputed:

  List returned by
  [`precompute_forward_inputs`](https://nickmckay.github.io/AARP/reference/precompute_forward_inputs.md).

- AAR_params:

  Named list with elements `Ae`, `Ea`, `R`, `x`.

- dT_yr:

  Integration timestep (yr). Default `1500`; coarse but accurate here
  because the abiotic accumulation is linear in time within each step
  (see package vignettes for a convergence check).

## Value

Numeric vector of predicted D/L, one per precomputed sample.
