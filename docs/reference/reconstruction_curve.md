# Reconstruction error curve

Reconstruction error curve

## Usage

``` r
reconstruction_curve(X, T, k = seq_len(ncol(T)), center = TRUE)
```

## Arguments

- X:

  A beat matrix.

- T:

  A basis matrix.

- k:

  Integer vector of numbers of components.

- center:

  Center before projecting.

## Value

A data frame with columns `k` and `error`.
