# Distance between two subspaces

Distance between two subspaces

## Usage

``` r
subspace_distance(T1, T2, type = c("rms", "max", "frobenius"))
```

## Arguments

- T1:

  A matrix whose columns span the first subspace.

- T2:

  A matrix whose columns span the second subspace.

- type:

  `"rms"` for the root mean squared sine of the principal angles,
  `"max"` for the largest angle, or `"frobenius"` for
  `||T1 T1' - T2 T2'||_F / sqrt(2 q)`.

## Value

A non-negative scalar.
