# PCA basis from an aggregated covariance matrix

Computes the eigenvectors of the average record covariance matrix. This
is the basis obtained by a federated PCA that shares `p x p` covariance
matrices instead of `p x q` bases.

## Usage

``` r
pool_covariance(
  X_list,
  q = 10L,
  weights = NULL,
  center = TRUE,
  return_cov = FALSE
)
```

## Arguments

- X_list:

  A list of beat matrices.

- q:

  Number of components.

- weights:

  Optional record weights.

- center:

  Center each record before computing its covariance.

- return_cov:

  Return the assembled covariance matrix as well.

## Value

A `bp_basis` object; with `return_cov = TRUE` a list with the basis `T`
and the covariance `S`.
