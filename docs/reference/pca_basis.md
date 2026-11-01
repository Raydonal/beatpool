# Principal component basis of a beat matrix

Principal component basis of a beat matrix

## Usage

``` r
pca_basis(X, q = 10L, center = TRUE, fix_sign = TRUE, scores = FALSE)
```

## Arguments

- X:

  A numeric matrix with beats in rows and time samples in columns.

- q:

  Number of components to keep.

- center:

  Center each time sample before the decomposition.

- fix_sign:

  Apply \[sign_fix()\] to the returned basis.

- scores:

  Return the projected scores as well.

## Value

A list with `T` (the `p x q` basis), `values` (explained variances),
`mean` (column means), `total` (total variance), `q` and, when
requested, `scores`.
