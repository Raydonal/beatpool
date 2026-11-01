# Simulate ECG-like heartbeat matrices

Generates one or several records of beats from a Gaussian model with a
low-rank covariance. Heterogeneity across records is introduced by
random rotations of a common eigenbasis.

## Usage

``` r
simulate_beats(
  H = 20L,
  n = 200L,
  p = 61L,
  q = 10L,
  model = c("homogeneous", "heterogeneous", "mixture"),
  epsilon = 0,
  mu = 1,
  pi = 0.5,
  sigma = 0.2,
  spectrum = NULL,
  basis = NULL,
  lambda_sd = 0,
  seed = NULL
)
```

## Arguments

- H:

  Number of records.

- n:

  Number of beats per record (a scalar or a vector of length `H`).

- p:

  Number of time samples.

- q:

  Number of signal components.

- model:

  `"homogeneous"`, `"heterogeneous"` or `"mixture"`.

- epsilon:

  Strength of the record-specific rotations.

- mu:

  Separation between the two class means in the mixture model.

- pi:

  Probability of the first class in the mixture model.

- sigma:

  Standard deviation of isotropic noise.

- spectrum:

  Optional vector of eigenvalues (length `q`).

- basis:

  Optional common `p x q` eigenbasis.

- lambda_sd:

  Standard deviation of log-normal record-specific eigenvalue
  perturbations.

- seed:

  Optional seed.

## Value

A list with the beat matrices `X`, the labels `labels`, the normal
indicators `normal`, the true record bases `Q_list`, the common basis
`Q0`, the eigenvalues `lambda`, and the simulation parameters.
