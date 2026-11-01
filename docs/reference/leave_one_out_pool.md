# Leave-one-record-out pooled bases

Leave-one-record-out pooled bases

## Usage

``` r
leave_one_out_pool(
  bases,
  data = NULL,
  method = "procrustes",
  align = TRUE,
  iterations = 1L,
  q = NULL
)
```

## Arguments

- bases:

  A list of bases.

- data:

  Optional list of beat matrices for column ordering.

- method:

  Pooling method.

- align:

  Align before pooling.

- iterations:

  Number of GPA iterations.

- q:

  Number of components.

## Value

A list with the full pooled basis `full` and a list `leave_one_out` of
pooled bases obtained by omitting each record.
