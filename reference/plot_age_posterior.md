# Plot the posterior age distribution for a two-acid age inversion

Plot the posterior age distribution for a two-acid age inversion

## Usage

``` r
plot_age_posterior(
  mcmc_out,
  burnin = 1000,
  true_age_ka = NULL,
  compare_ages_ka = NULL
)
```

## Arguments

- mcmc_out:

  List returned by
  [`invert_two_acid_age`](https://nickmckay.github.io/AARP/reference/invert_two_acid_age.md).

- burnin:

  Integer. Burn-in iterations to discard. Default `1000`.

- true_age_ka:

  Optional numeric. Draws a reference vertical line (e.g. an
  independently known age, for validation). Default `NULL`.

- compare_ages_ka:

  Optional named numeric vector of point age estimates to compare
  against (e.g. single-acid or GLS-combined ages), drawn as dashed
  vertical lines. Default `NULL`.

## Value

A ggplot object.
