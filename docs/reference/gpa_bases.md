# Generalised Procrustes analysis of bases

Iteratively aligns and averages a list of bases until the pooled
subspace stabilises.

## Usage

``` r
gpa_bases(
  bases,
  method = "procrustes",
  q = NULL,
  data = NULL,
  weights = NULL,
  tol = 1e-10,
  max_iter = 50L
)
```

## Arguments

- bases:

  A list of matrices or `bp_basis` objects.

- method:

  One of `"procrustes"`, `"bar"`, `"qr"` or `"mean"`.

- q:

  Number of components to pool. Defaults to the smallest number of
  columns among the inputs.

- data:

  Optional list of beat matrices used to order the pooled columns by
  average explained variance.

- weights:

  Optional record weights.

- tol:

  Convergence tolerance for the iterative algorithm.

- max_iter:

  Maximum number of iterations.

## Value

A `bp_pool` object.
