# Compress and reconstruct beats

Compress and reconstruct beats

## Usage

``` r
compress_beats(X, T, k = ncol(T), center = TRUE)
```

## Arguments

- X:

  A beat matrix.

- T:

  A basis matrix.

- k:

  Number of components.

- center:

  Center before projecting.

## Value

A list with the coefficients `Y`, the reconstructed beats `X_hat`, the
per-beat relative errors `error` and the mean relative error
`mean_error`.
