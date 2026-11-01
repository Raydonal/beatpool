# Fix the sign of the columns of a basis

Multiplies each column by -1 when its largest absolute entry is
negative. This makes the basis reproducible without changing its span.

## Usage

``` r
sign_fix(T)
```

## Arguments

- T:

  A matrix with orthonormal (or arbitrary) columns.

## Value

The sign-adjusted matrix.
