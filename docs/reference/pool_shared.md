# Several shared bases from transmitted moments

Partitions the records into `C` groups and returns one basis per group.
The receiver alternates between the leading `k` eigenvectors of each
group sum and the assignment that retains the most of each moment, from
`starts` random assignments. It uses only what the records have already
transmitted, so no further communication from them is needed.

## Usage

``` r
pool_shared(moments, k, C = 3L, starts = 10L, iters = 50L, seed = NULL)
```

## Arguments

- moments:

  A list of `bp_moment` objects.

- k:

  Dimension of each basis.

- C:

  Number of shared bases.

- starts:

  Number of random initial assignments.

- iters:

  Maximum number of alternations per start.

- seed:

  Optional seed. The random state of the session is restored.

## Value

An object of class `bp_shared` with the list of bases `T`, the
assignment `group` and the objective `value`.

## Examples

``` r
sim <- simulate_beats(H = 6, n = 50, p = 20, q = 3, seed = 3)
m <- lapply(sim$X, energy_moment, q = 5)
fit <- pool_shared(m, k = 2, C = 2, starts = 2, seed = 1)
fit$group
#> [1] 2 2 1 2 2 1
```
