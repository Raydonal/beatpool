# Reproducing the article

These scripts reproduce every table, figure and number quoted in the article
on shared heartbeat compression bases built from transmitted second moments.
They use only exported functions of `beatpool`, base R, `ggplot2` and `MASS`.

```r
Sys.setenv(BP_DATA = "data", BP_OUT = "article-output", BP_CORES = "4")
source(system.file("article", "run_all.R", package = "beatpool"))
```

| Step | Script | Output in `BP_OUT` |
|---|---|---|
| 1 | `01_data.R` | downloads the MIT-BIH Arrhythmia, MIT-BIH Supraventricular Arrhythmia and St Petersburg INCART databases from PhysioNet, the 49 industry portfolios of the Kenneth R. French data library and the FRED recession indicator; writes `beats_*.rds` and `corpora.csv` |
| 2 | `02_heldout.R` | held-out errors, rate--distortion grid, basis choice and the normal-versus-ventricular task (`<corpus>_*.csv`) |
| 3 | `03_summary.R` | paired bootstrap intervals (`rules.csv`), means, rates, LaTeX tables and `values.tex` |
| 4 | `04_figures.R` | the figures of the article (PDF) |

`values.tex` defines, through `\bpval{key}{text}`, every number quoted in the
text of the article, so that the manuscript reads them instead of copying them.

The PhysioNet databases are fixed releases (version 1.0.0). The French data
library revises its returns from time to time, so a fresh download may change
the financial results in the third decimal place; the article used the file
downloaded in October 2026.
