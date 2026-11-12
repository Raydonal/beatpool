# Choose attention settings from transmitted moments

Leave-one-out over the records, using only what they have transmitted:
each record in turn is the query, attends to the others, and the
retained error \\1 - tr(T^T M_j T)\\ of its own rank-q moment is
averaged over records. Because \\M_j\\ is truncated, this error is an
upper bound on the per-beat error of the record. By default the query of
record \\j\\ is \\M_j\\ itself, which is less noisy than a query
computed from a few beats and favours smaller temperatures; `queries`
gives queries computed as they will be for a new record, such as moments
of the first beats, at the cost of one more transmitted moment per
record.

## Usage

``` r
select_attention(moments, k, grid, exclude = NULL, queries = NULL)
```

## Arguments

- moments:

  A list of `bp_moment` objects.

- k:

  Dimension of the keys.

- grid:

  A data frame with one row per candidate and columns among `tau`,
  `score`, `normalize`, `sink`, `contrast`, `self`, `top`, `heads` (a
  comma-separated string) and `refine` (the temperature of
  \[refine_moments()\]). `NA` in `sink`, `top`, `heads` or `refine`
  leaves that option off. Refinement is recomputed without the left-out
  record.

- exclude:

  Optional list giving, for each record, the indices of the records it
  must not attend to, such as other records of the same patient.

- queries:

  Optional list of query matrices, one per record.

## Value

A list with the chosen row `best` and the mean error of every row,
`errors`.

## Examples

``` r
sim <- simulate_beats(H = 6, n = 60, p = 20, q = 3, seed = 5)
m <- lapply(sim$X, energy_moment, q = 5)
g <- expand.grid(tau = c(0.01, 0.1, Inf), sink = c(NA, 0))
select_attention(m, k = 3, grid = g)$best
#>    tau sink
#> 1 0.01   NA
```
