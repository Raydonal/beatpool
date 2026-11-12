## Step 1. Download the four corpora and build the beat matrices.
## PhysioNet databases are read in WFDB format; the financial panel comes from
## the Kenneth R. French data library and the FRED recession indicator.

source(file.path(Sys.getenv("BP_ARTICLE", system.file("article", package = "beatpool")), "helpers.R"))

## ECG databases: first signal of MIT-BIH and SVDB, lead II of INCART.
ecg <- list(
  mitbih = list(db = "mitdb", dir = "mitdb", wsize = 30L, lead = NULL),
  svdb = list(db = "svdb", dir = "svdb", wsize = 15L, lead = NULL),
  incart = list(db = "incartdb", dir = "incart", wsize = 21L, lead = "II"))
for (t in names(ecg)) {
  e <- ecg[[t]]
  d <- file.path(data_dir, e$dir)
  download_physionet(e$db, d)
  X <- read_database(d, wsize = e$wsize, lead = e$lead, min_beats = 100L)
  saveRDS(X, out(sprintf("beats_%s.rds", t)))
  message(CORPORA[[t]], ": ", length(X), " records")
}

## Financial panel: 49 value-weighted industry portfolios, daily, 1970 to 2025,
## one group per calendar year with at least 200 days, standardised in the year.
fin <- file.path(data_dir, "french")
dir.create(fin, showWarnings = FALSE, recursive = TRUE)
csv <- file.path(fin, "49_Industry_Portfolios_Daily.csv")
if (!file.exists(csv)) {
  zip <- file.path(fin, "49_Industry_Portfolios_daily_CSV.zip")
  utils::download.file(paste0("https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/",
                              "ftp/49_Industry_Portfolios_daily_CSV.zip"), zip, mode = "wb")
  utils::unzip(zip, exdir = fin)
}
lines <- readLines(csv, warn = FALSE)
i0 <- grep("Average Value Weighted Returns", lines)[1L]
industries <- strsplit(lines[i0 + 1L], ",")[[1]][-1]
i1 <- i0 + 2L
while (i1 <= length(lines) && grepl("^[0-9]{8}", lines[i1])) i1 <- i1 + 1L
dat <- read.csv(text = lines[(i0 + 2L):(i1 - 1L)], header = FALSE)
names(dat) <- c("date", industries)
dat$date <- as.Date(as.character(dat$date), format = "%Y%m%d")
R <- as.matrix(dat[, industries])
R[R <= -99] <- NA_real_
ok <- stats::complete.cases(R) & dat$date >= as.Date("1970-01-01") &
  dat$date <= as.Date("2025-12-31")
dat <- dat[ok, ]
R <- R[ok, , drop = FALSE]

fred <- file.path(data_dir, "fred", "USREC.csv")
if (!file.exists(fred)) {
  dir.create(dirname(fred), showWarnings = FALSE, recursive = TRUE)
  utils::download.file("https://fred.stlouisfed.org/graph/fredgraph.csv?id=USREC", fred)
}
rec <- read.csv(fred)
names(rec) <- c("date", "usrec")
recession <- rec$usrec[match(format(dat$date, "%Y-%m"), format(as.Date(rec$date), "%Y-%m"))]
recession[is.na(recession)] <- 0L
label <- ifelse(recession == 1L, "recession", "expansion")

year <- format(dat$date, "%Y")
panel <- list()
for (y in sort(unique(year))) {
  i <- which(year == y)
  if (length(i) < 200L) next
  Z <- scale(R[i, , drop = FALSE])
  attr(Z, "record") <- y
  attr(Z, "labels") <- label[i]
  attr(Z, "normal") <- label[i] == "expansion"
  panel[[y]] <- Z
}
saveRDS(panel, out("beats_finance.rds"))

## Table 1: groups, beats, dimension and the mean share of variance outside
## ten components.
corpus_table <- do.call(rbind, lapply(names(CORPORA), function(t) {
  X <- readRDS(out(sprintf("beats_%s.rds", t)))
  tail10 <- vapply(X, function(x) {
    b <- pca_basis(x, q = 10L)
    1 - sum(b$values) / b$total
  }, numeric(1))
  groups <- if (t == "incart") length(unique(vapply(X, attr, "", "patient"))) else length(X)
  data.frame(corpus = CORPORA[[t]], records = length(X), groups = groups,
             beats = sum(vapply(X, nrow, integer(1))), p = ncol(X[[1]]),
             tail10 = mean(tail10))
}))
write.csv(corpus_table, out("corpora.csv"), row.names = FALSE)
print(corpus_table)
