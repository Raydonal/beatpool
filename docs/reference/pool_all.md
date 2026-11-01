# All pooling operators at once

Runs the operator families compared in the paper on a common set of
local bases and returns their \\p \times q\\ operators, together with
the number of values each operator transmits per record. The function is
the programmatic counterpart of the comparison tables: the returned
matrices can be passed to \[reconstruction_error()\] or
\[subspace_distance()\].

## Usage

``` r
pool_all(
  X_list = NULL,
  bases = NULL,
  q = 10L,
  weights = NULL,
  residual = FALSE
)
```

## Arguments

- X_list:

  A list of beat matrices, one per record. Required unless `bases` is
  supplied.

- bases:

  Optional list of `bp_basis` objects. When `NULL`, bases are extracted
  from `X_list`.

- q:

  Number of components.

- weights:

  Optional record weights.

- residual:

  Use residual imputation for the spectral summaries.

## Value

A list of class `bp_poolall`. Each element is a \\p \times q\\ matrix
with orthonormal columns, except the unaligned means, which are
column-normalised but not orthogonal. The attribute `numbers` gives the
communication cost of each operator per record.

## Examples

``` r
sim <- simulate_beats(H = 6, n = 120, p = 20, q = 3, seed = 1)
all <- pool_all(sim$X, q = 3)
names(all)
#>  [1] "bar"  "bara" "cov"  "covn" "gpa"  "hier" "proc" "qr"   "qra"  "spec"
attr(all, "numbers")
#>  bar bara  cov covn  gpa hier proc   qr  qra spec 
#>   60   60  210  210   60   60   60   60   60   63 
```
