# Changelog

## beatpool (development version)

- [`simulate_beats()`](https://raydonal.github.io/beatpool/reference/simulate_beats.md)
  labels every beat `"normal"` (with `normal = TRUE`) in the homogeneous
  and heterogeneous models, which have a single class; they were `NA`.
  The simulated signals are unchanged.

- Raydonal Ospina is now the sole author and copyright holder; Abel
  Borges is acknowledged in the README.

- The package now reproduces every table, figure and number of the
  accompanying article: the scripts in
  `system.file("article", package = "beatpool")` download the data and
  run the analysis (see the vignette “Reproducing the article”).

- New
  [`download_physionet()`](https://raydonal.github.io/beatpool/reference/download_physionet.md)
  downloads a PhysioNet database, and
  [`read_database()`](https://raydonal.github.io/beatpool/reference/read_database.md)
  reads every record of a WFDB database into beat matrices with the
  settings of the article.

- New attention over transmitted moments.
  [`attention_weights()`](https://raydonal.github.io/beatpool/reference/attention_weights.md)
  and
  [`pool_attention()`](https://raydonal.github.io/beatpool/reference/pool_attention.md)
  score each transmitted moment against a query moment of the beats to
  be compressed and pool the moments with softmax weights; equal weights
  are the limit of infinite temperature. Options: symmetric projection
  and cosine scores, sparsemax weights, a sink key holding the mean
  moment, a subtracted mean moment (differential attention),
  interpolation with the query moment, several heads and sparse weights.

- New
  [`select_attention()`](https://raydonal.github.io/beatpool/reference/select_attention.md)
  and
  [`select_tau()`](https://raydonal.github.io/beatpool/reference/select_tau.md):
  settings chosen by leave-one-out over the transmitted moments,
  optionally with queries of the size that a new record will have.

- New
  [`refine_moments()`](https://raydonal.github.io/beatpool/reference/refine_moments.md):
  self-attention among records with a symmetric, doubly stochastic
  matrix found by Sinkhorn iterations; it leaves the mean moment
  unchanged.

- New
  [`stream_attention()`](https://raydonal.github.io/beatpool/reference/stream_attention.md):
  one basis per block of beats, with the moment of the earlier beats as
  query.

## beatpool 0.2.0

- New
  [`energy_moment()`](https://raydonal.github.io/beatpool/reference/energy_moment.md)
  and
  [`pool_energy()`](https://raydonal.github.io/beatpool/reference/pool_energy.md):
  moments of beats normalised by their energy, optimal for the mean
  per-beat relative error.
- New
  [`pool_shared()`](https://raydonal.github.io/beatpool/reference/pool_shared.md)
  and
  [`select_shared()`](https://raydonal.github.io/beatpool/reference/select_shared.md):
  several shared bases fitted from the moments already transmitted by
  the records.

## beatpool 0.1.0

- First version: moment pooling, direction pooling, hierarchical
  pooling, residual imputation, metrics, simulation and a WFDB reader in
  R.
