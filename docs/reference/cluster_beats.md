# Cluster beats in a projected basis

Cluster beats in a projected basis

## Usage

``` r
cluster_beats(
  X,
  T,
  k = 3L,
  K = 2L,
  truth = NULL,
  center = TRUE,
  nstart = 10L,
  iter.max = 100L,
  seed = NULL
)
```

## Arguments

- X:

  A beat matrix.

- T:

  A basis matrix.

- k:

  Number of leading components used as features.

- K:

  Number of clusters.

- truth:

  Optional true labels used to compute quality metrics.

- center:

  Center before projecting.

- nstart:

  Number of random starts.

- iter.max:

  Maximum number of k-means iterations.

- seed:

  Optional seed.

## Value

A list with the cluster labels, the fitted `kmeans` object, the features
`Y` and, when `truth` is supplied, `metrics`.
