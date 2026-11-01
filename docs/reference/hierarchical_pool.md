# Hierarchical Procrustes pooling

Pools the leading `k`-dimensional subspaces separately for every
`k = 1, ..., q`. Each subspace is aligned to a medoid reference before
averaging, which avoids forcing a single compromise basis to serve all
dimensions.

## Usage

``` r
hierarchical_pool(
  bases,
  q = 10L,
  ref = "medoid",
  iterations = 1L,
  data = NULL,
  robust = FALSE,
  ...
)
```

## Arguments

- bases:

  A list of matrices or `bp_basis` objects.

- q:

  Maximum number of components.

- ref:

  Reference rule: `"medoid"` (default) or an integer index. The medoid
  is computed once in the full `q`-dimensional subspace and the same
  reference record is used for every dimension.

- iterations:

  Number of GPA iterations for each subspace.

- data:

  Optional list of beat matrices used to order components.

- robust:

  Use robust reweighting.

- ...:

  Passed to \[pool_bases()\].

## Value

A list of class `bp_hpool` with one `bp_pool` object per dimension.
