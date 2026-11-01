# Read a WFDB annotation file

Reads MIT-BIH style annotation files. Each annotation is one 16-bit word
containing a 10-bit time increment and a 6-bit code. Auxiliary strings
(code 63) and skip records (code 59) are handled.

## Usage

``` r
read_annotations(file)
```

## Arguments

- file:

  Path to a `.atr` (or other annotation) file.

## Value

A data frame with one row per annotation and columns `sample` (0-based
WFDB sample), `index` (1-based sample), `code`, `symbol`, `beat` and
`normal`.
