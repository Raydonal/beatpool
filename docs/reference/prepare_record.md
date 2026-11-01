# Prepare a record from a signal and its annotations

Prepare a record from a signal and its annotations

## Usage

``` r
prepare_record(
  x,
  ann,
  record = "record",
  fs = 360,
  wsize = 30L,
  lead = "MLII",
  normalize = c("segment", "qrs", "none"),
  drop_first_last = TRUE
)
```

## Arguments

- x:

  Numeric ECG signal.

- ann:

  Annotation data frame produced by \[read_annotations()\].

- record:

  Record identifier.

- fs:

  Sampling frequency.

- wsize:

  Half-window size.

- lead:

  Lead name.

- normalize:

  Normalisation method passed to \[segment_beats()\].

- drop_first_last:

  Drop the first and last annotated beats.

## Value

A matrix of beats with attributes, as returned by \[segment_beats()\].
