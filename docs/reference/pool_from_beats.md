# Pooled basis from a list of beat matrices

Convenience wrapper that extracts the bases and pools them in one call.

## Usage

``` r
pool_from_beats(
  X_list,
  q = 10L,
  method = "procrustes",
  align = TRUE,
  iterations = 1L,
  weights = NULL,
  ...
)
```

## Arguments

- X_list:

  A list of beat matrices.

- q:

  Number of components.

- method:

  Pooling method.

- align:

  Align before pooling.

- iterations:

  Number of GPA iterations.

- weights:

  Optional record weights.

- ...:

  Passed to \[pool_bases()\].

## Value

A `bp_pool` object.
