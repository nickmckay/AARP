# Simulate a placeholder glacial-interglacial temperature history

The calibration dataset records only a single "representative" bottom
water temperature per sample – not an actual paleotemperature
reconstruction. This generates a mild synthetic oscillation around that
recorded value, giving the depth-series time integration a temperature
*history* to work with rather than a constant. The amplitude and period
are placeholders (order-of-magnitude plausible for glacial-interglacial
deep polar ocean variability), not a real reconstruction, and are meant
to be replaced with actual per-core paleotemperature histories when
available.

## Usage

``` r
simulate_temp_history(
  temp_C,
  max_age_ka,
  amplitude = 1,
  period_ka = 100,
  resolution_ka = 5
)
```

## Arguments

- temp_C:

  Recorded (present-day / representative) bottom water temperature (deg
  C).

- max_age_ka:

  Oldest age (ka BP) the history needs to cover.

- amplitude:

  Oscillation amplitude (deg C). Default `1`.

- period_ka:

  Oscillation period (ka). Default `100` (approximate Pleistocene
  glacial-interglacial pacing).

- resolution_ka:

  Spacing of control points (ka). Default `5`.

## Value

A data frame with columns `age_ka`, `temp_C`, suitable for
[`get_temp_at_time`](https://nickmckay.github.io/AARP/reference/get_temp_at_time.md)
/
[`racemize_one_depth`](https://nickmckay.github.io/AARP/reference/racemize_one_depth.md).
