# Package index

## Bases

Basis extraction and construction.

- [`pca_basis()`](https://raydonal.github.io/beatpool/reference/pca_basis.md)
  : Principal component basis of a beat matrix

- [`dhmt_basis()`](https://raydonal.github.io/beatpool/reference/dhmt_basis.md)
  : Discrete Hermite transform basis

- [`dct_basis()`](https://raydonal.github.io/beatpool/reference/dct_basis.md)
  : Discrete cosine transform basis (DCT-II)

- [`random_basis()`](https://raydonal.github.io/beatpool/reference/random_basis.md)
  : Random orthogonal basis

- [`extract_bases()`](https://raydonal.github.io/beatpool/reference/extract_bases.md)
  : Extract PCA bases from a list of beat matrices

- [`sign_fix()`](https://raydonal.github.io/beatpool/reference/sign_fix.md)
  : Fix the sign of the columns of a basis

- [`orthonormalize()`](https://raydonal.github.io/beatpool/reference/orthonormalize.md)
  :

  Orthonormalise the first `k` columns of a basis

- [`complete_basis()`](https://raydonal.github.io/beatpool/reference/complete_basis.md)
  : Complete a partial basis to a full orthogonal matrix

- [`order_basis()`](https://raydonal.github.io/beatpool/reference/order_basis.md)
  : Order the columns of a basis by explained variance

- [`basis_variance()`](https://raydonal.github.io/beatpool/reference/basis_variance.md)
  : Average explained variance of each basis column

- [`basis_frame()`](https://raydonal.github.io/beatpool/reference/basis_frame.md)
  : Convert a basis to a tidy data frame

## Alignment and geometry

Removing the indeterminacy of principal components.

- [`align_basis()`](https://raydonal.github.io/beatpool/reference/align_basis.md)
  : Align a basis to a reference basis
- [`principal_angles()`](https://raydonal.github.io/beatpool/reference/principal_angles.md)
  : Principal angles between two subspaces
- [`subspace_distance()`](https://raydonal.github.io/beatpool/reference/subspace_distance.md)
  : Distance between two subspaces

## Pooling

Combining bases across records and sites.

- [`pool_bases()`](https://raydonal.github.io/beatpool/reference/pool_bases.md)
  : Pool patient-specific bases
- [`gpa_bases()`](https://raydonal.github.io/beatpool/reference/gpa_bases.md)
  : Generalised Procrustes analysis of bases
- [`hierarchical_pool()`](https://raydonal.github.io/beatpool/reference/hierarchical_pool.md)
  : Hierarchical Procrustes pooling
- [`pool_spectral()`](https://raydonal.github.io/beatpool/reference/pool_spectral.md)
  : Pool second moments from per-record spectral summaries
- [`pool_covariance()`](https://raydonal.github.io/beatpool/reference/pool_covariance.md)
  : PCA basis from an aggregated covariance matrix
- [`energy_moment()`](https://raydonal.github.io/beatpool/reference/energy_moment.md)
  : Energy-normalised moment of one record
- [`pool_energy()`](https://raydonal.github.io/beatpool/reference/pool_energy.md)
  : Pool energy-normalised moments
- [`pool_shared()`](https://raydonal.github.io/beatpool/reference/pool_shared.md)
  : Several shared bases from transmitted moments
- [`select_shared()`](https://raydonal.github.io/beatpool/reference/select_shared.md)
  : Choose a shared basis for a record
- [`attention_weights()`](https://raydonal.github.io/beatpool/reference/attention_weights.md)
  : Attention weights over transmitted moments
- [`pool_attention()`](https://raydonal.github.io/beatpool/reference/pool_attention.md)
  : Attention pooling of transmitted moments
- [`select_tau()`](https://raydonal.github.io/beatpool/reference/select_tau.md)
  : Choose the attention temperature from transmitted moments
- [`select_attention()`](https://raydonal.github.io/beatpool/reference/select_attention.md)
  : Choose attention settings from transmitted moments
- [`refine_moments()`](https://raydonal.github.io/beatpool/reference/refine_moments.md)
  : Self-attention among records with doubly stochastic weights
- [`stream_attention()`](https://raydonal.github.io/beatpool/reference/stream_attention.md)
  : Attention with a query that grows with the record
- [`pool_from_beats()`](https://raydonal.github.io/beatpool/reference/pool_from_beats.md)
  : Pooled basis from a list of beat matrices
- [`pool_all()`](https://raydonal.github.io/beatpool/reference/pool_all.md)
  : All pooling operators at once
- [`leave_one_out_pool()`](https://raydonal.github.io/beatpool/reference/leave_one_out_pool.md)
  : Leave-one-record-out pooled bases

## Reconstruction and metrics

- [`reconstruction_error()`](https://raydonal.github.io/beatpool/reference/reconstruction_error.md)
  : Relative reconstruction error
- [`reconstruction_curve()`](https://raydonal.github.io/beatpool/reference/reconstruction_curve.md)
  : Reconstruction error curve
- [`compress_beats()`](https://raydonal.github.io/beatpool/reference/compress_beats.md)
  : Compress and reconstruct beats
- [`rate_distortion()`](https://raydonal.github.io/beatpool/reference/rate_distortion.md)
  : Rate–distortion curve of a quantised spectral codec
- [`ari()`](https://raydonal.github.io/beatpool/reference/ari.md) :
  Adjusted Rand index
- [`purity()`](https://raydonal.github.io/beatpool/reference/purity.md)
  : Cluster purity
- [`prop_correct()`](https://raydonal.github.io/beatpool/reference/prop_correct.md)
  : Local cluster consistency
- [`cluster_metrics()`](https://raydonal.github.io/beatpool/reference/cluster_metrics.md)
  : Cluster quality metrics
- [`cluster_beats()`](https://raydonal.github.io/beatpool/reference/cluster_beats.md)
  : Cluster beats in a projected basis
- [`jackknife_stability()`](https://raydonal.github.io/beatpool/reference/jackknife_stability.md)
  : Jackknife stability of a pooled basis
- [`bootstrap_ci()`](https://raydonal.github.io/beatpool/reference/bootstrap_ci.md)
  : Percentile bootstrap confidence interval
- [`paired_bootstrap()`](https://raydonal.github.io/beatpool/reference/paired_bootstrap.md)
  : Paired bootstrap comparison
- [`compare_methods()`](https://raydonal.github.io/beatpool/reference/compare_methods.md)
  : Compare several methods across records

## Simulation

- [`simulate_beats()`](https://raydonal.github.io/beatpool/reference/simulate_beats.md)
  : Simulate ECG-like heartbeat matrices
- [`corrupt_basis()`](https://raydonal.github.io/beatpool/reference/corrupt_basis.md)
  : Corrupt a basis with random signs and permutations

## ECG preprocessing and WFDB

- [`segment_beats()`](https://raydonal.github.io/beatpool/reference/segment_beats.md)
  : Segment ECG heartbeats around R peaks

- [`prepare_record()`](https://raydonal.github.io/beatpool/reference/prepare_record.md)
  : Prepare a record from a signal and its annotations

- [`read_mitdb_csv()`](https://raydonal.github.io/beatpool/reference/read_mitdb_csv.md)
  :

  Read a MIT-BIH CSV exported by `rdsamp`

- [`read_mitdb_record()`](https://raydonal.github.io/beatpool/reference/read_mitdb_record.md)
  : Read and prepare a MIT-BIH record from CSV and annotation files

- [`read_wfdb_header()`](https://raydonal.github.io/beatpool/reference/read_wfdb_header.md)
  : Read a WFDB header file

- [`read_wfdb()`](https://raydonal.github.io/beatpool/reference/read_wfdb.md)
  : Read a WFDB signal record

- [`download_physionet()`](https://raydonal.github.io/beatpool/reference/download_physionet.md)
  : Download a PhysioNet database

- [`read_database()`](https://raydonal.github.io/beatpool/reference/read_database.md)
  : Beat matrices of a WFDB database

- [`read_annotations()`](https://raydonal.github.io/beatpool/reference/read_annotations.md)
  : Read a WFDB annotation file

- [`ann_codes`](https://raydonal.github.io/beatpool/reference/ann_codes.md)
  : WFDB annotation codes

- [`beat_summary()`](https://raydonal.github.io/beatpool/reference/beat_summary.md)
  : Beat-level summaries

## Plotting

- [`plot_beats()`](https://raydonal.github.io/beatpool/reference/plot_beats.md)
  : Plot a sample of heartbeats
- [`plot_mean_beat()`](https://raydonal.github.io/beatpool/reference/plot_mean_beat.md)
  : Plot the mean heartbeat
- [`plot_basis()`](https://raydonal.github.io/beatpool/reference/plot_basis.md)
  : Plot basis coefficients
- [`plot_scree()`](https://raydonal.github.io/beatpool/reference/plot_scree.md)
  : Plot a scree or cumulative variance curve
- [`plot_reconstruction()`](https://raydonal.github.io/beatpool/reference/plot_reconstruction.md)
  : Plot reconstruction error curves
- [`plot_cost_error()`](https://raydonal.github.io/beatpool/reference/plot_cost_error.md)
  : Plot communication cost against reconstruction error
- [`plot_clustering()`](https://raydonal.github.io/beatpool/reference/plot_clustering.md)
  : Plot clustering accuracy curves
- [`plot_stability()`](https://raydonal.github.io/beatpool/reference/plot_stability.md)
  : Plot jackknife stability
- [`plot_simulation()`](https://raydonal.github.io/beatpool/reference/plot_simulation.md)
  : Plot simulation results
- [`theme_bp()`](https://raydonal.github.io/beatpool/reference/theme_bp.md)
  : Consistent ggplot2 theme
- [`method_palette()`](https://raydonal.github.io/beatpool/reference/method_palette.md)
  : Colour palette for pooling methods
