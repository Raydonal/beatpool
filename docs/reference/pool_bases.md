# Pool patient-specific bases

Combines a list of orthonormal bases into a single basis. Four pooling
operators are available: the arithmetic mean, the column-normalised
mean, the QR retraction and the orthogonal Procrustes mean. The
Procrustes mean is the minimiser of the sum of squared Frobenius
distances to the inputs and is the recommended default.

## Usage

``` r
pool_bases(
  bases,
  method = c("procrustes", "bar", "qr", "mean"),
  align = TRUE,
  ref = 1L,
  iterations = 1L,
  q = NULL,
  data = NULL,
  weights = NULL,
  tol = 1e-10,
  max_iter = 50L,
  robust = FALSE,
  robust_c = 1.345
)
```

## Arguments

- bases:

  A list of matrices or `bp_basis` objects.

- method:

  One of `"procrustes"`, `"bar"`, `"qr"` or `"mean"`.

- align:

  Align the bases before pooling. Alignment removes the sign,
  permutation and rotation indeterminacy of principal components.

- ref:

  Reference used by the one-step alignment: an integer index into
  `bases`, or a matrix. Ignored when `iterations > 1`.

- iterations:

  Number of alignment/pooling iterations. Values greater than one run
  generalised Procrustes analysis.

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

- robust:

  Apply robust reweighting to downweight outlying records.

- robust_c:

  Tuning constant for the robust weights.

## Value

A list of class `bp_pool` with the pooled basis `T`, the method, the
aligned bases, the per-record rotations, the average explained variances
and convergence information.
