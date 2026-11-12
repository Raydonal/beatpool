# Attention weights over transmitted moments

Scores each transmitted moment against a query moment and turns the
scores into weights. The key of record \\h\\ is the projector \\P_h\\
onto its leading `k` eigenvectors. With `score = "energy"` the score is
\\s_h = tr(P_h R_0)\\, where \\R_0\\ is the energy-normalised moment of
the query; when \\R_0\\ is the full moment of a set of beats, \\1 -
s_h\\ is their mean per-beat relative error in the basis of record
\\h\\. With `score = "projection"` the score is \\tr(P_h P_0)/k\\,
symmetric in the two subspaces, where \\P_0\\ projects onto the leading
`k` eigenvectors of \\R_0\\. With `score = "cosine"` it is the cosine
between \\R_0\\ and \\M_h\\ in the Frobenius inner product. The scores
depend on the eigenvectors only through \\P_h\\ or \\M_h\\, and are
therefore unchanged by their signs, order or any rotation within the
retained subspace.

## Usage

``` r
attention_weights(
  query,
  moments,
  k,
  tau,
  top = NULL,
  score = "energy",
  normalize = "softmax",
  sink = NULL
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

## Value

A vector of weights over the records, with the scores in the attribute
`scores`, the weight of the extra key in `sink` and the effective number
of records \\1 / \sum_h w_h^2\\ in `n_eff`.

## Examples

``` r
sim <- simulate_beats(H = 5, n = 60, p = 20, q = 3, seed = 4)
m <- lapply(sim$X, energy_moment, q = 5)
attention_weights(sim$X[[1]][1:20, ], m[-1], k = 3, tau = 0.01)
#> [1] 0.1374985 0.1566100 0.1343125 0.5715790
#> attr(,"scores")
#> [1] 0.6172103 0.6185117 0.6169759 0.6314582
#> attr(,"sink")
#> [1] 0
#> attr(,"n_eff")
#> [1] 2.576158
```
