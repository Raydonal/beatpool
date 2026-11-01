# Pool second moments from per-record spectral summaries

Each record is summarised by its leading eigenpairs \\(Q_h, \lambda_h)\\
and the pooled covariance is assembled as \$\$\hat\Sigma = \sum_h w_h
Q_h \Lambda_h Q_h'.\$\$ The outer product is invariant to signs,
permutations and rotations of \\Q_h\\, so no alignment is needed. When
`residual = TRUE` the variance outside the retained subspace is spread
over its orthogonal complement, which keeps the trace of each record
covariance intact.

## Usage

``` r
pool_spectral(
  bases,
  q = 10L,
  q_out = q,
  weights = NULL,
  residual = FALSE,
  return_cov = FALSE
)
```

## Arguments

- bases:

  A list of `bp_basis` objects.

- q:

  Number of leading eigenpairs used in the per-record summaries.

- q_out:

  Number of components returned. It may exceed `q`.

- weights:

  Optional record weights.

- residual:

  Impute the residual variance over the discarded subspace.

- return_cov:

  Return the assembled covariance matrix as well.

## Value

A `bp_basis` object; with `return_cov = TRUE` a list with the basis `T`
and the covariance `S`.
