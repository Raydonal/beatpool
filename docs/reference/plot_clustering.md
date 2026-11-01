# Plot clustering accuracy curves

Plot clustering accuracy curves

## Usage

``` r
plot_clustering(df, metric = "Adjusted Rand index")
```

## Arguments

- df:

  A data frame with columns `k`, `value`, `method` and optionally `K`.

- metric:

  Label for the vertical axis.

## Value

A ggplot object.
