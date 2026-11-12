# Beat matrices of a WFDB database

Reads every record of a database stored in WFDB format, segments one
channel around the annotated beats and keeps records with at least
`min_beats` beats. These are the settings of the accompanying article:
the first signal with `wsize = 30` for the MIT-BIH Arrhythmia Database
and `wsize = 15` for the MIT-BIH Supraventricular Arrhythmia Database,
and lead II with `wsize = 21` for the St Petersburg INCART database.

## Usage

``` r
read_database(dir, wsize, lead = NULL, min_beats = 100L, normalize = "segment")
```

## Arguments

- dir:

  Directory with the `.hea`, `.dat` and `.atr` files.

- wsize:

  Half window in samples; a beat has `2 * wsize + 1` samples.

- lead:

  Signal description to read, or `NULL` for the first signal.

- min_beats:

  Minimum number of beats of a record.

- normalize:

  Passed to \[prepare_record()\].

## Value

A named list of beat matrices. For records whose header names a patient,
the attribute `patient` holds it.

## Examples

``` r
if (FALSE) { # \dontrun{
download_physionet("svdb", "data/svdb")
X <- read_database("data/svdb", wsize = 15)
} # }
```
