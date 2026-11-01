# Relative reconstruction error

Projects each beat onto the subspace spanned by the first `k` columns of
a basis and returns the relative squared error.

## Usage

``` r
reconstruction_error(X, T, k = ncol(T), center = TRUE, per_beat = FALSE)
```

## Arguments

- X:

  A beat matrix.

- T:

  A basis matrix.

- k:

  Number of leading components.

- center:

  Center each column of `X` before projecting.

- per_beat:

  Return one value per beat instead of the average.

## Value

A numeric scalar or vector.
