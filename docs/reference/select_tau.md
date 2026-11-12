# Choose the attention temperature from transmitted moments

A shortcut for \[select_attention()\] over temperatures alone.

## Usage

``` r
select_tau(
  moments,
  k,
  taus = c(0, 2^seq(-10, 0), Inf),
  exclude = NULL,
  queries = NULL
)
```

## Arguments

- moments:

  A list of `bp_moment` objects.

- k:

  Dimension of the keys.

- taus:

  Candidate temperatures, `Inf` included by default.

- exclude:

  Optional list giving, for each record, the indices of the records it
  must not attend to, such as other records of the same patient.

- queries:

  Optional list of query matrices, one per record.

## Value

The chosen temperature, with the mean error of each candidate in the
attribute `errors`.

## Examples

``` r
sim <- simulate_beats(H = 6, n = 60, p = 20, q = 3, seed = 5)
m <- lapply(sim$X, energy_moment, q = 5)
select_tau(m, k = 3)
#> [1] 0.0078125
#> attr(,"errors")
#>            0 0.0009765625  0.001953125   0.00390625    0.0078125     0.015625 
#>    0.3388116    0.3345384    0.3322226    0.3302791    0.3297223    0.3298872 
#>      0.03125       0.0625        0.125         0.25          0.5            1 
#>    0.3301638    0.3303607    0.3304743    0.3305349    0.3305662    0.3305820 
#>          Inf 
#>    0.3305980 
```
