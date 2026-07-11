# Build per-sample age-depth models from grouped core data

Groups calibration rows by `(study, core)` and constructs an age-depth
model for each row. Cores with two or more distinct depths get a real
interpolated
[`make_age_model`](https://nickmckay.github.io/AARP/reference/make_age_model.md)
(duplicate depths are collapsed by averaging their ages first); cores
with only one usable depth fall back to
[`rate_to_age_model`](https://nickmckay.github.io/AARP/reference/rate_to_age_model.md),
a straight line from the core top (0 cm, 0 ka) through that single
point.

## Usage

``` r
build_core_age_models(calibration_data)
```

## Arguments

- calibration_data:

  Data frame with columns `study`, `core`, `depth_mbsf`, `age_ka` (one
  row per calibration sample).

## Value

A list, one element per row of `calibration_data` (same order), each an
age-depth model as returned by
[`make_age_model`](https://nickmckay.github.io/AARP/reference/make_age_model.md)
/
[`rate_to_age_model`](https://nickmckay.github.io/AARP/reference/rate_to_age_model.md).
