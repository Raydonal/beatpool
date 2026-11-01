# Choose a shared basis for a record

Returns the index of the basis with the lowest mean per-beat relative
error on the given beats. A compressor transmits this index, of
`log2(C)` bits.

## Usage

``` r
select_shared(fit, X)
```

## Arguments

- fit:

  A `bp_shared` object.

- X:

  A beat matrix.

## Value

An integer index, with the errors in the attribute `errors`.

## Examples

``` r
sim <- simulate_beats(H = 6, n = 50, p = 20, q = 3, seed = 3)
m <- lapply(sim$X, energy_moment, q = 5)
fit <- pool_shared(m, k = 2, C = 2, starts = 2, seed = 1)
select_shared(fit, sim$X[[1]])
#> [1] 2
#> attr(,"errors")
#> [1] 0.5584512 0.5101670
```
