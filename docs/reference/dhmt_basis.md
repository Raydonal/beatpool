# Discrete Hermite transform basis

Builds the eigenvectors of the tridiagonal matrix used by Mugler's
discrete Hermite transform. The construction follows the one used in the
original ECG analysis.

## Usage

``` r
dhmt_basis(p, b = 1)
```

## Arguments

- p:

  Number of time samples.

- b:

  Bandwidth parameter.

## Value

A `p x p` orthogonal matrix.
