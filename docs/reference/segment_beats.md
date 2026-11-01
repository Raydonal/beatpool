# Segment ECG heartbeats around R peaks

Extracts fixed-length windows around annotated R peaks. Each window has
`2 * wsize + 1` samples. By default the window is standardised with the
mean and standard deviation of the whole beat segment delimited by the
midpoints between consecutive R peaks, which is the preprocessing used
in the original MIT-BIH analysis.

## Usage

``` r
segment_beats(
  x,
  peaks,
  labels = NULL,
  normal = NULL,
  record = "record",
  fs = 360,
  wsize = 30L,
  normalize = c("segment", "qrs", "none"),
  drop_edges = TRUE,
  drop_first_last = FALSE
)
```

## Arguments

- x:

  Numeric ECG signal.

- peaks:

  Integer vector of R-peak positions (1-based).

- labels:

  Optional beat labels, one per peak.

- normal:

  Optional logical vector, one per peak.

- record:

  Record identifier stored in the result.

- fs:

  Sampling frequency in Hz.

- wsize:

  Half-window size. The window length is `2 * wsize + 1`.

- normalize:

  One of `"segment"` (default), `"qrs"` or `"none"`.

- drop_edges:

  Drop beats whose window falls outside the signal.

- drop_first_last:

  Drop the first and last annotated beats.

## Value

A matrix of beats (rows) by time samples (columns), with attributes
`peaks`, `labels`, `normal`, `record`, `fs`, `wsize` and `normalize`.
