# Percentile bootstrap confidence interval

Percentile bootstrap confidence interval

## Usage

``` r
bootstrap_ci(x, B = 2000L, conf = 0.95, statistic = mean, seed = NULL)
```

## Arguments

- x:

  Numeric vector.

- B:

  Number of bootstrap replicates.

- conf:

  Confidence level.

- statistic:

  Function returning the statistic of interest.

- seed:

  Optional seed.

## Value

A named numeric vector with the lower limit, estimate and upper limit.
