# Orthonormalise the first `k` columns of a basis

Returns an orthonormal basis for the span of the selected columns. This
puts non-orthogonal pooling rules such as the column-normalised mean on
the same footing as orthogonal ones when computing projections or
clustering features.

## Usage

``` r
orthonormalize(T, k = ncol(T), fix_sign = TRUE)
```

## Arguments

- T:

  A basis matrix.

- k:

  Number of leading columns to use.

- fix_sign:

  Apply \[sign_fix()\] to the result.

## Value

A `p x k` matrix with orthonormal columns.
