## Reproduce every table, figure and number of the article.
##
##   Rscript -e 'source(system.file("article", "run_all.R", package = "beatpool"))'
##
## Environment variables: BP_DATA (downloads, default "data"), BP_OUT (results,
## default "article-output") and BP_CORES (default 4). With 4 cores the whole
## run takes about 20 minutes after the downloads. Set OPENBLAS_NUM_THREADS=1
## when the BLAS is multithreaded.

dir <- Sys.getenv("BP_ARTICLE", system.file("article", package = "beatpool"))
Sys.setenv(BP_ARTICLE = dir)
for (s in c("01_data.R", "02_heldout.R", "03_summary.R", "04_figures.R")) {
  message("== ", s)
  source(file.path(dir, s), local = new.env())
}
