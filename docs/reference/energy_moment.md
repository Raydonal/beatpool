# Energy-normalised moment of one record

Centres the beats, divides each one by its norm and returns the leading
eigenpairs of the resulting moment \\R_h = n_h^{-1} \sum_i x\_{hi}
x\_{hi}^T / \\x\_{hi}\\^2\\. For a matrix \\T\\ with orthonormal
columns, the mean per-beat relative error of the record is \\1 - tr(T^T
R_h T)\\, so that the moment targets that error directly. The sample
covariance targets the ratio of traces instead.

## Usage

``` r
energy_moment(X, q = 10L)
```

## Arguments

- X:

  A beat matrix with beats in rows.

- q:

  Number of eigenpairs retained.

## Value

An object of class `bp_moment` with the eigenvectors `T`, the
eigenvalues `values`, the unit-trace tail `tail` and the number of beats
`n`.

## Examples

``` r
sim <- simulate_beats(H = 2, n = 50, p = 20, q = 3, seed = 1)
m <- energy_moment(sim$X[[1]], q = 5)
m$tail
#> [1] 0.2802386
```
