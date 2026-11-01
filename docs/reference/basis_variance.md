# Average explained variance of each basis column

Average explained variance of each basis column

## Usage

``` r
basis_variance(T, X_list, weights = NULL)
```

## Arguments

- T:

  A basis matrix.

- X_list:

  A list of beat matrices used to score the columns.

- weights:

  Optional record weights.

## Value

A numeric vector with one value per column of `T`.
