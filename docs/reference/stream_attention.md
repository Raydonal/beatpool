# Attention with a query that grows with the record

Splits the beats into consecutive blocks. The basis of each block after
the first comes from attention with the energy-normalised moment of all
earlier beats as query, so that it uses no beat of the block itself. The
transmitted moments play the part of a fixed cache of keys and values,
and only the query is updated.

## Usage

``` r
stream_attention(X, moments, k, block = 50L, args = list(tau = Inf))
```

## Arguments

- X:

  A beat matrix, in time order.

- moments:

  A list of `bp_moment` objects.

- k:

  Dimension of the bases.

- block:

  Number of beats per block.

- args:

  Arguments of \[pool_attention()\] other than the query, the moments
  and `k`, as a list, or a function of the number of earlier beats
  returning such a list.

## Value

A list with one basis per block, `T`, and the block of each beat,
`block`. The first block has no basis.

## Examples

``` r
sim <- simulate_beats(H = 5, n = 60, p = 20, q = 3, seed = 7)
m <- lapply(sim$X[-1], energy_moment, q = 5)
s <- stream_attention(sim$X[[1]], m, k = 3, block = 20,
                      args = list(tau = 0.05))
length(s$T)
#> [1] 3
```
