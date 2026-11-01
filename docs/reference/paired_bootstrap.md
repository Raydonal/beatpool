# Paired bootstrap comparison

Paired bootstrap comparison

## Usage

``` r
paired_bootstrap(x, y, B = 2000L, conf = 0.95, seed = NULL)
```

## Arguments

- x:

  Numeric vector for method one.

- y:

  Numeric vector for method two, paired with `x`.

- B:

  Number of bootstrap replicates.

- conf:

  Confidence level.

- seed:

  Optional seed.

## Value

A list with the observed mean difference, its confidence interval and a
two-sided bootstrap p-value.
