# Read a WFDB header file

Parses the `.hea` file of a WFDB record. Formats 212, 16 and 80 are
supported by \[read_wfdb()\].

## Usage

``` r
read_wfdb_header(file)
```

## Arguments

- file:

  Path to a `.hea` file.

## Value

A list with the record name, sampling frequency, number of samples, and
a data frame describing each signal.
