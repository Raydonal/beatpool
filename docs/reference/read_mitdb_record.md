# Read and prepare a MIT-BIH record from CSV and annotation files

Read and prepare a MIT-BIH record from CSV and annotation files

## Usage

``` r
read_mitdb_record(
  record,
  data_dir = ".",
  lead = c("MLII", "V5"),
  wsize = 30L,
  normalize = c("segment", "qrs", "none"),
  drop_first_last = TRUE
)
```

## Arguments

- record:

  Record name, e.g. `"100"`.

- data_dir:

  Directory containing `<record>.csv` and `atr/<record>.atr`.

- lead:

  ECG lead to use.

- wsize:

  Half-window size.

- normalize:

  Normalisation method.

- drop_first_last:

  Drop the first and last annotated beats.

## Value

A matrix of beats with attributes.
