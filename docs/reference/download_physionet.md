# Download a PhysioNet database

Downloads the header, signal and annotation files of every record listed
in the `RECORDS` file of a PhysioNet database. Files already present in
`dir` are kept.

## Usage

``` r
download_physionet(
  db,
  dir,
  version = "1.0.0",
  ext = c("hea", "dat", "atr"),
  records = NULL
)
```

## Arguments

- db:

  Database slug, such as `"mitdb"`, `"svdb"` or `"incartdb"`.

- dir:

  Destination directory.

- version:

  Database version.

- ext:

  File extensions downloaded for each record.

- records:

  Optional subset of records.

## Value

The record names, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
download_physionet("mitdb", "data/mitdb")
} # }
```
