# Align a basis to a reference basis

Finds the orthogonal matrix `O` that maximises
`trace(t(ref) %*% Q %*% O)`. This removes the sign, permutation and
rotation indeterminacy of principal component bases before they are
pooled.

## Usage

``` r
align_basis(Q, ref)
```

## Arguments

- Q:

  A basis to align.

- ref:

  The reference basis.

## Value

A list with the aligned basis `Q`, the rotation `O` and the achieved
correlation `correlation`.
