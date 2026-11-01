# Order the columns of a basis by explained variance

Order the columns of a basis by explained variance

## Usage

``` r
order_basis(T, X_list, weights = NULL, fix_sign = TRUE)
```

## Arguments

- T:

  A basis matrix.

- X_list:

  A list of beat matrices.

- weights:

  Optional record weights.

- fix_sign:

  Apply \[sign_fix()\] after reordering.

## Value

A list with the reordered basis `T` and the variances `values`.
