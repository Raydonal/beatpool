# Self-attention among records with doubly stochastic weights

Replaces each transmitted moment by an average of all of them, \\\tilde
M_j = \sum_h A\_{jh} M_h\\, with \\A = D K D\\, \\K\_{jh} = \exp\\tr(P_j
P_h)/(k\tau)\\\\ and \\D\\ diagonal, found by Sinkhorn iterations so
that every row and column of \\A\\ sums to one. Because \\K\\ is
symmetric, so is \\A\\: each record gives as much weight as it receives,
and no record dominates the refined set. The receiver can do this with
what it holds.

## Usage

``` r
refine_moments(moments, k, tau, iters = 500L)
```

## Arguments

- moments:

  A list of `bp_moment` objects.

- k:

  Dimension of the projectors.

- tau:

  Temperature.

- iters:

  Maximum number of Sinkhorn iterations.

## Value

A list of `bp_moment` objects of the same rank, with the matrix `A` in
the attribute `A`.

## Examples

``` r
sim <- simulate_beats(H = 5, n = 60, p = 20, q = 3, seed = 6)
m <- lapply(sim$X, energy_moment, q = 5)
r <- refine_moments(m, k = 3, tau = 0.05)
round(colSums(attr(r, "A")), 6)
#> [1] 1 1 1 1 1
```
