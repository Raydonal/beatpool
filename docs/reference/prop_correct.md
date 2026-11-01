# Local cluster consistency

For each beat, the proportion of correct assignments among the other
beats in the same cluster must exceed one half. This is the criterion
used in the original ECG analysis.

## Usage

``` r
prop_correct(cluster, truth)
```

## Arguments

- cluster:

  Cluster labels.

- truth:

  True labels.

## Value

The mean proportion of locally consistent beats.
