# beatpool

**Moment pooling of PCA bases for ECG heartbeats**

`beatpool` represents heartbeats with patient-specific principal
component bases, pools those bases across sites without sharing raw
signals, and provides the reconstruction, clustering and stability tools
used in the accompanying study. A pure-R reader for WFDB headers,
signals and annotations is included, so the full pipeline runs without
external WFDB software.

## Why pool bases?

Each recording has its own covariance structure, so a basis fitted on
one patient is not directly comparable with a basis fitted on another.
Pooling the bases entry by entry fails because principal components are
defined up to signs, permutations and rotations. `beatpool` offers two
families of solutions:

- **Alignment-based pooling** removes the indeterminacy with an
  orthogonal Procrustes rotation and then averages, or takes the
  orthogonal Procrustes mean (the default), the column-normalised mean
  or a QR retraction.
- **Moment pooling** never aligns anything. Each site sends its leading
  eigenpairs \\(Q_h, \lambda_h)\\ and the server assembles \\\hat\Sigma
  = \sum_h w_h Q_h \Lambda_h Q_h'\\. The outer product is invariant to
  any orthogonal transformation of \\Q_h\\, and the leading eigenvectors
  of \\\hat\Sigma\\ are Eckart-Young optimal for the pooled second
  moment.

## Installation

``` r
# from a local clone
install.packages("remotes")
remotes::install_local("beatpool", dependencies = TRUE)
```

## Quick start

``` r
library(beatpool)

sim <- simulate_beats(H = 8, n = 200, p = 61, q = 10,
                      model = "heterogeneous", epsilon = 0.2, seed = 1)
bases <- extract_bases(sim$X, q = 10)

pooled <- pool_bases(bases, method = "procrustes")   # aligned mean
spec   <- pool_spectral(bases, q = 10)               # moment pooling

reconstruction_error(sim$X[[1]], pooled$T, k = 10)
reconstruction_error(sim$X[[1]], spec$T, k = 10)

# every operator in the paper, with its communication cost
all <- pool_all(sim$X, q = 10)
attr(all, "numbers")
```

Reading a MIT-BIH record requires only the `.hea`, `.dat` and `.atr`
files:

``` r
x <- read_wfdb("100", dir = "data/mitdb")
ann <- read_annotations("data/mitdb/atr/100.atr")
X <- prepare_record(x$MLII, ann, record = "100", fs = 360, wsize = 30)
b <- pca_basis(X, q = 10)
```

## Main functions

| Task           | Functions                                                                                                                                                                                                                                                                                                                                                                                                                           |
|----------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Bases          | [`pca_basis()`](https://raydonal.github.io/beatpool/reference/pca_basis.md), [`dhmt_basis()`](https://raydonal.github.io/beatpool/reference/dhmt_basis.md), [`dct_basis()`](https://raydonal.github.io/beatpool/reference/dct_basis.md), [`random_basis()`](https://raydonal.github.io/beatpool/reference/random_basis.md)                                                                                                          |
| Alignment      | [`align_basis()`](https://raydonal.github.io/beatpool/reference/align_basis.md), [`principal_angles()`](https://raydonal.github.io/beatpool/reference/principal_angles.md), [`subspace_distance()`](https://raydonal.github.io/beatpool/reference/subspace_distance.md)                                                                                                                                                             |
| Pooling        | [`pool_bases()`](https://raydonal.github.io/beatpool/reference/pool_bases.md), [`pool_spectral()`](https://raydonal.github.io/beatpool/reference/pool_spectral.md), [`hierarchical_pool()`](https://raydonal.github.io/beatpool/reference/hierarchical_pool.md), [`pool_covariance()`](https://raydonal.github.io/beatpool/reference/pool_covariance.md), [`pool_all()`](https://raydonal.github.io/beatpool/reference/pool_all.md) |
| Reconstruction | [`reconstruction_error()`](https://raydonal.github.io/beatpool/reference/reconstruction_error.md), [`compress_beats()`](https://raydonal.github.io/beatpool/reference/compress_beats.md), [`plot_cost_error()`](https://raydonal.github.io/beatpool/reference/plot_cost_error.md)                                                                                                                                                   |
| Clustering     | [`cluster_beats()`](https://raydonal.github.io/beatpool/reference/cluster_beats.md), [`ari()`](https://raydonal.github.io/beatpool/reference/ari.md), [`purity()`](https://raydonal.github.io/beatpool/reference/purity.md), [`prop_correct()`](https://raydonal.github.io/beatpool/reference/prop_correct.md)                                                                                                                      |
| Stability      | [`leave_one_out_pool()`](https://raydonal.github.io/beatpool/reference/leave_one_out_pool.md), [`jackknife_stability()`](https://raydonal.github.io/beatpool/reference/jackknife_stability.md), [`bootstrap_ci()`](https://raydonal.github.io/beatpool/reference/bootstrap_ci.md)                                                                                                                                                   |
| Simulation     | [`simulate_beats()`](https://raydonal.github.io/beatpool/reference/simulate_beats.md), [`corrupt_basis()`](https://raydonal.github.io/beatpool/reference/corrupt_basis.md)                                                                                                                                                                                                                                                          |
| WFDB           | [`read_wfdb()`](https://raydonal.github.io/beatpool/reference/read_wfdb.md), [`read_annotations()`](https://raydonal.github.io/beatpool/reference/read_annotations.md), [`prepare_record()`](https://raydonal.github.io/beatpool/reference/prepare_record.md)                                                                                                                                                                       |

## Documentation

- [`vignette("beatpool")`](https://raydonal.github.io/beatpool/articles/beatpool.md)
  walks through simulation, pooling, reconstruction and clustering.
- The pkgdown site documents every exported function.
- The manuscript in `paper/` derives the invariance and optimality
  results and reports the MIT-BIH, SVDB and financial experiments.

## Citation

``` r
citation("beatpool")
```

## Acknowledgements

The author thanks Abel Borges, whose undergraduate presentation on ECG
compression with averaged PCA bases, given more than eight years ago,
motivated this work.

## License

MIT (c) Raydonal Ospina.
