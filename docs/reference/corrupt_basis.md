# Corrupt a basis with random signs and permutations

Used to demonstrate the need for alignment before pooling.

## Usage

``` r
corrupt_basis(Q, mode = c("sign", "permutation", "both"), seed = NULL)
```

## Arguments

- Q:

  A basis matrix.

- mode:

  One of `"sign"`, `"permutation"` or `"both"`.

- seed:

  Optional seed.

## Value

A corrupted basis with the same span.
