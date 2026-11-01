# Extract PCA bases from a list of beat matrices

Extract PCA bases from a list of beat matrices

## Usage

``` r
extract_bases(X_list, q = 10L, ...)
```

## Arguments

- X_list:

  A list of beat matrices, one per record.

- q:

  Number of components to keep.

- ...:

  Passed to \[pca_basis()\].

## Value

A list of `bp_basis` objects.
