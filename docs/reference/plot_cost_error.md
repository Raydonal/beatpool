# Plot communication cost against reconstruction error

Places a set of operators in the cost–error plane. The cost is read from
the `numbers` attribute of a `bp_poolall` object, or supplied directly,
and the error is a named vector indexed by operator.

## Usage

``` r
plot_cost_error(x, error, labels = names(error), k = 3L)
```

## Arguments

- x:

  A `bp_poolall` object or a named numeric vector of costs.

- error:

  A named numeric vector of reconstruction errors.

- labels:

  Optional labels, one per operator. Defaults to the names.

- k:

  Number of coefficients, used in the axis label only.

## Value

A ggplot object.

## Examples

``` r
sim <- simulate_beats(H = 6, n = 120, p = 20, q = 3, seed = 2)
all <- pool_all(sim$X, q = 3)
err <- vapply(all, function(T) mean(vapply(sim$X, reconstruction_error,
  numeric(1), T = T, k = 3)), numeric(1))
plot_cost_error(all, err)
```
