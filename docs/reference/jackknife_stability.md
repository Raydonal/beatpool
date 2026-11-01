# Jackknife stability of a pooled basis

Jackknife stability of a pooled basis

## Usage

``` r
jackknife_stability(loo, type = "rms")
```

## Arguments

- loo:

  A list returned by \[leave_one_out_pool()\].

- type:

  Distance used between the full and leave-one-out bases.

## Value

A data frame with the relative Frobenius distance and the subspace
distance for each omitted record.
