# Read a MIT-BIH CSV exported by `rdsamp`

Read a MIT-BIH CSV exported by `rdsamp`

## Usage

``` r
read_mitdb_csv(file, lead = c("MLII", "V5"))
```

## Arguments

- file:

  Path to the `.csv` file.

- lead:

  One of `"MLII"` or `"V5"`.

## Value

A data frame with `time` and `x`.
