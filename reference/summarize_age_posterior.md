# Posterior summary for a two-acid age inversion

Posterior summary for a two-acid age inversion

## Usage

``` r
summarize_age_posterior(mcmc_out, burnin = 1000)
```

## Arguments

- mcmc_out:

  List returned by
  [`invert_two_acid_age`](https://nickmckay.github.io/AARP/reference/invert_two_acid_age.md).

- burnin:

  Integer. Burn-in iterations to discard. Default `1000`.

## Value

A one-row data frame with columns `median_ka`, `lo50_ka`, `hi50_ka`,
`lo66_ka`, `hi66_ka`, `lo95_ka`, `hi95_ka`.
