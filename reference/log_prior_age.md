# Log-prior for sample age

A weak, bounded-uniform prior on \\\log(\text{age})\\, wide enough to
cover the full calibration range (default `c(0.5, 3000)` ka) without
favouring any age within it. Because this is an ordinary Bayesian prior,
it can be tightened using independent information (e.g. a stratigraphic
age range, or a known minimum/maximum age) simply by narrowing
`range_ka` or swapping in an informative density – something the
original GLS point-estimate calculator cannot easily accommodate.

## Usage

``` r
log_prior_age(log_age, range_ka = c(0.5, 3000))
```

## Arguments

- log_age:

  Candidate value of \\\log(\text{age in ka})\\.

- range_ka:

  Numeric length-2 vector, the allowed age range (ka). Default
  `c(0.5, 3000)`.

## Value

Scalar log-prior (`0` inside the range, `-Inf` outside).
