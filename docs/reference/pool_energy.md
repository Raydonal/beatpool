# Pool energy-normalised moments

Diagonalises the average of the rank-q moments transmitted by the
records. The cost per record is the same as that of \[pool_spectral()\].

## Usage

``` r
pool_energy(moments, q_out = moments[[1]]$q, weights = NULL)
```

## Arguments

- moments:

  A list of `bp_moment` objects.

- q_out:

  Number of components returned.

- weights:

  Optional record weights.

## Value

A `bp_basis` object.

## Examples

``` r
sim <- simulate_beats(H = 4, n = 50, p = 20, q = 3, seed = 2)
m <- lapply(sim$X, energy_moment, q = 5)
pool_energy(m, q_out = 3)$T[1:2, ]
#>            [,1]        [,2]        [,3]
#> [1,] 0.01281598 -0.02407280 0.002527693
#> [2,] 0.01081224  0.01535701 0.002659275
```
