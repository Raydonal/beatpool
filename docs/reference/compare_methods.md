# Compare several methods across records

Compare several methods across records

## Usage

``` r
compare_methods(values, reference = colnames(values)[1])
```

## Arguments

- values:

  A data frame with one row per record and one column per method.

- reference:

  Name of the reference method.

## Value

A data frame with the mean, bootstrap interval and paired comparison
against the reference for each method.
