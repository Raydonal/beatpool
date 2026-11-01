# Read a WFDB signal record

Pure-R reader for WFDB signal files in formats 212, 16 and 80. No
external WFDB software is required.

## Usage

``` r
read_wfdb(record, dir = ".", channels = NULL, from = 1L, to = NULL)
```

## Arguments

- record:

  Record name (e.g. `"100"`).

- dir:

  Directory containing the `.hea` and `.dat` files.

- channels:

  Integer vector of channels to return; `NULL` returns all.

- from:

  First sample (1-based) to return.

- to:

  Last sample (1-based) to return; `NULL` returns to the end.

## Value

A data frame with a `time` column (seconds) and one column per channel,
named after the header descriptions. The parsed header is stored in the
`"header"` attribute.
