# Attention pooling of transmitted moments

Returns the leading eigenvectors of an attended moment \\S = \sum_h w_h
M_h\\, where \\M_h\\ are the rank-q moments already transmitted and
\\w_h\\ are the \[attention_weights()\] of a query. With `tau = Inf` and
the defaults it coincides with \[pool_energy()\]; with `tau = 0` it
returns the basis of the record of highest score. No communication from
the records is added; the query is computed where the beats to be
compressed are.

## Usage

``` r
pool_attention(
  query,
  moments,
  k,
  tau,
  q_out = k,
  top = NULL,
  score = "energy",
  normalize = "softmax",
  sink = NULL,
  contrast = 0,
  self = 0,
  heads = NULL
)
```

## Arguments

- query:

  A `bp_moment` object, a \\p \times p\\ moment matrix or a beat matrix,
  from which a full-rank energy moment is computed.

- moments:

  A list of `bp_moment` objects.

- k:

  Dimension of the keys.

- tau:

  Temperature. `Inf` gives equal weights and `0` puts all weight on the
  highest score, shared among ties.

- q_out:

  Number of components returned. Defaults to `k`.

- top:

  Optional number of largest weights kept, the others set to zero and
  the kept ones renormalised. A compressor then transmits `top` indices
  and weights instead of one weight per record.

- score:

  One of `"energy"`, `"projection"` or `"cosine"`.

- normalize:

  `"softmax"`, or `"sparsemax"`, the Euclidean projection of \\s/\tau\\
  onto the simplex, which sets low scores to zero.

- sink:

  Optional offset of an extra key, the mean of the transmitted moments,
  whose score is its own score plus `sink`. Its weight, in the attribute
  `sink`, goes to the population moment.

- contrast:

  Weight \\\lambda\\ of the subtracted mean moment.

- self:

  Weight \\\alpha\\ of the query moment.

- heads:

  Optional vector of key dimensions, one per head.

## Value

A list with the basis `T`, the eigenvalues `values`, the `weights` over
the records (a matrix with one row per head when there are several), the
weight of the extra key `sink`, the temperature `tau` and `n_eff`.

## Details

Three options change the attended moment, each term scaled to unit trace
first. `contrast` subtracts \\\lambda\\ times the mean moment, as in
differential attention, so that the directions shared by all records
weigh less. `self` mixes in the query, \\(1-\alpha) S + \alpha R_0\\, as
in retrieval interpolation; the basis then depends on the record's own
beats and must be sent with them. `heads` averages the attended moments
of several key dimensions.

## Examples

``` r
sim <- simulate_beats(H = 5, n = 60, p = 20, q = 3, seed = 4)
m <- lapply(sim$X, energy_moment, q = 5)
fit <- pool_attention(sim$X[[1]][1:20, ], m[-1], k = 3, tau = 0.01)
round(fit$weights, 3)
#> [1] 0.137 0.157 0.134 0.572
```
