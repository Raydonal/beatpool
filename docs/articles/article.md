# Reproducing the article

The package ships the scripts that reproduce every table, figure and
number of the accompanying article. They download the public data, so
they are not run when the vignette is built.

``` r
Sys.setenv(BP_DATA = "data", BP_OUT = "article-output", BP_CORES = "4",
           OPENBLAS_NUM_THREADS = "1")
source(system.file("article", "run_all.R", package = "beatpool"))
```

## What each step does

1.  `01_data.R` calls
    [`download_physionet()`](https://raydonal.github.io/beatpool/reference/download_physionet.md)
    for the MIT-BIH Arrhythmia Database, the MIT-BIH Supraventricular
    Arrhythmia Database and the St Petersburg INCART database, and
    [`read_database()`](https://raydonal.github.io/beatpool/reference/read_database.md)
    to segment the first signal (lead II in INCART) into windows of 61,
    31 and 43 samples around the annotated R peaks. It also downloads
    the 49 industry portfolios of the French data library and
    standardizes them within each year.
2.  `02_heldout.R` leaves out one record at a time, or one patient in
    INCART, and estimates every operator from the others:
    [`pool_spectral()`](https://raydonal.github.io/beatpool/reference/pool_spectral.md),
    [`pool_covariance()`](https://raydonal.github.io/beatpool/reference/pool_covariance.md),
    [`pool_bases()`](https://raydonal.github.io/beatpool/reference/pool_bases.md),
    [`hierarchical_pool()`](https://raydonal.github.io/beatpool/reference/hierarchical_pool.md),
    [`energy_moment()`](https://raydonal.github.io/beatpool/reference/energy_moment.md)
    with
    [`pool_shared()`](https://raydonal.github.io/beatpool/reference/pool_shared.md),
    and
    [`dct_basis()`](https://raydonal.github.io/beatpool/reference/dct_basis.md).
    It records the mean per-beat error, the ratio of traces, the rate of
    [`rate_distortion()`](https://raydonal.github.io/beatpool/reference/rate_distortion.md)
    and the separation of normal and ventricular beats.
3.  `03_summary.R` computes paired bootstrap intervals over records or
    patients and writes the tables and `values.tex`.
4.  `04_figures.R` draws the figures.

The scripts live in `system.file("article", package = "beatpool")`.
