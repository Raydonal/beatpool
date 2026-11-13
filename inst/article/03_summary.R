## Step 3. Paired bootstrap intervals, mean errors, rates and the downstream
## task. Writes CSV summaries, the LaTeX tables of the article and
## values.tex, from which the article takes every number quoted in the text.

source(file.path(Sys.getenv("BP_ARTICLE", system.file("article", package = "beatpool")), "helpers.R"))
B <- 2000L
set.seed(14)
tags <- names(CORPORA)
rd <- function(t, f) read.csv(out(sprintf("%s_%s.csv", t, f)), colClasses = c(record = "character"))
groups <- lapply(tags, function(t) {
  d <- read.csv(out(sprintf("%s_clusters.csv", t)), colClasses = c(record = "character", group = "character"))
  setNames(d$group, d$record)
})
names(groups) <- tags

## Paired mean difference a - b over groups, with a bootstrap over groups.
paired <- function(a, b, g) {
  d <- tapply(a - b, g, mean)
  bs <- replicate(B, mean(sample(d, length(d), replace = TRUE)))
  c(delta = mean(d), lo = unname(quantile(bs, 0.025)), hi = unname(quantile(bs, 0.975)),
    median = unname(median(d)))
}
wide <- function(h, val) {
  w <- reshape(h[, c("record", "method", val)], idvar = "record", timevar = "method",
               direction = "wide")
  names(w) <- sub(paste0("^", val, "\\."), "", names(w))
  w
}
pairs <- list(c("multi3", "spec"), c("enorm", "spec"), c("multi2", "spec"), c("multi4", "spec"),
              c("multi3_n0", "spec"), c("spec", "covn"), c("spec", "proc"), c("spec", "hier"),
              c("multi3", "enorm"), c("multi3", "proc"), c("multi3", "hier"), c("enorm", "hier"),
              c("multi3", "local"))
rules <- do.call(rbind, lapply(tags, function(t) {
  h <- rd(t, "heldout")
  do.call(rbind, lapply(c(3L, 10L), function(k) do.call(rbind, lapply(c("beat", "trace"), function(v) {
    w <- wide(h[h$k == k, ], v)
    do.call(rbind, lapply(pairs, function(z) {
      x <- paired(w[[z[1]]], w[[z[2]]], groups[[t]][w$record])
      data.frame(corpus = t, k = k, metric = v, a = z[1], b = z[2], t(x),
                 better = mean(w[[z[1]]] < w[[z[2]]]))
    }))
  }))))
}))
write.csv(rules, out("rules.csv"), row.names = FALSE)

means <- do.call(rbind, lapply(tags, function(t) {
  h <- rd(t, "heldout")
  h$g <- groups[[t]][h$record]
  a <- aggregate(cbind(beat, trace) ~ method + k + g, h, mean)
  a <- aggregate(cbind(beat, trace) ~ method + k, a, mean)
  a$corpus <- t
  a
}))
write.csv(means, out("means.csv"), row.names = FALSE)

## Bits per beat at given distortions. Shared bases cost 16 bits per entry,
## charged once over all beats of the corpus (broadcast) or to each record.
TARGETS <- c(20, 15, 10, 7, 5)
inv <- function(d, target, col) {
  d <- d[order(d[[col]]), ]
  y <- cummin(d$D)
  if (target < min(y) || target > max(y)) return(NA_real_)
  stats::approx(y, d[[col]], xout = target, rule = 1, ties = mean)$y
}
p_of <- function(t) ncol(readRDS(out(sprintf("beats_%s.rds", t)))[[1]])
rates <- do.call(rbind, lapply(tags, function(t) {
  r <- rd(t, "rate")
  tot <- sum(unique(r[, c("record", "n")])$n)
  bits <- 10 * p_of(t) * 16
  r$index <- ifelse(r$method == "multi3", log2(3) / r$n, 0)
  r$per_record <- r$rate0 + r$nbasis * bits / r$n + r$index
  r$broadcast <- r$rate0 + ifelse(r$method == "local", bits / r$n, r$nbasis * bits / tot) + r$index
  a <- aggregate(cbind(D, per_record, broadcast) ~ method + step, r, mean)
  write.csv(a, out(sprintf("%s_rate_curve.csv", t)), row.names = FALSE)
  do.call(rbind, lapply(unique(a$method), function(m) {
    d <- a[a$method == m, ]
    do.call(rbind, lapply(TARGETS, function(x) data.frame(
      corpus = t, method = m, target = x, floor = min(d$D),
      per_record = inv(d, x, "per_record"), broadcast = inv(d, x, "broadcast"))))
  }))
}))
write.csv(rates, out("rates.csv"), row.names = FALSE)

down <- do.call(rbind, lapply(c("mitbih", "svdb", "incart"), function(t) {
  a <- aggregate(cbind(bal, lda_bal) ~ record + method + k, rd(t, "downstream"), mean)
  do.call(rbind, lapply(c(3L, 10L), function(k) do.call(rbind, lapply(c("bal", "lda_bal"), function(v) {
    w <- wide(a[a$k == k, ], v)
    do.call(rbind, lapply(c("multi3", "enorm", "local", "proc"), function(m) {
      x <- paired(w[[m]], w$spec, groups[[t]][w$record])
      data.frame(corpus = t, k = k, metric = v, method = m, n = nrow(w),
                 spec = mean(w$spec), value = mean(w[[m]]), t(x))
    }))
  }))))
}))
write.csv(down, out("downstream.csv"), row.names = FALSE)

chosen <- do.call(rbind, lapply(tags, function(t) {
  d <- rd(t, "chosen")
  data.frame(corpus = t, k = d$k, agree = d$chosen == d$chosen_n0)
}))
selection <- aggregate(agree ~ corpus + k, chosen, mean)
write.csv(selection, out("selection.csv"), row.names = FALSE)

## values.tex: \bpval{key}{text}
f3 <- function(x) gsub("-", "$-$", sprintf("%.3f", x))
f3m <- function(x) sprintf("$%s$", sub("^-", "-", sprintf("%.3f", x)))
z3 <- function(x) { x <- round(x, 3); x[x == 0] <- 0; x }  # no negative zero
ci <- function(r) sprintf("$%.3f$ [$%.3f$, $%.3f$]", z3(r$delta), z3(r$lo), z3(r$hi))
pct <- function(x) sprintf("%.0f\\%%", 100 * x)
val <- character()
put <- function(key, text) val[[key]] <<- text
for (i in seq_len(nrow(rules))) {
  r <- rules[i, ]
  key <- paste(r$corpus, r$a, r$b, r$k, r$metric, sep = ":")
  put(paste0("d:", key), ci(r))
  put(paste0("delta:", key), sprintf("$%.3f$", z3(r$delta)))
  put(paste0("median:", key), sprintf("$%.3f$", z3(r$median)))
  put(paste0("better:", key), pct(r$better))
  put(paste0("worse:", key), pct(1 - r$better))
}
for (i in seq_len(nrow(means))) {
  r <- means[i, ]
  put(paste("m", r$corpus, r$method, r$k, sep = ":"), sprintf("$%.3f$", r$beat))
  put(paste("mt", r$corpus, r$method, r$k, sep = ":"), sprintf("$%.3f$", r$trace))
}
for (i in seq_len(nrow(rates))) {
  r <- rates[i, ]
  put(paste("floor", r$corpus, r$method, sep = ":"), sprintf("%.1f\\%%", r$floor))
  if (!is.na(r$broadcast)) put(paste("bits", r$corpus, r$method, r$target, sep = ":"), sprintf("$%.1f$", r$broadcast))
  if (!is.na(r$per_record)) put(paste("bitsrec", r$corpus, r$method, r$target, sep = ":"), sprintf("$%.1f$", r$per_record))
}
for (i in seq_len(nrow(down))) {
  r <- down[i, ]
  put(paste("ds", r$corpus, r$metric, r$method, r$k, sep = ":"), ci(r))
}
for (i in seq_len(nrow(selection)))
  put(paste("agree", selection$corpus[i], selection$k[i], sep = ":"), pct(selection$agree[i]))
cp <- read.csv(out("corpora.csv"))
cp$tag <- tags
q <- 10
for (i in seq_len(nrow(cp))) {
  t <- cp$tag[i]; p <- cp$p[i]
  put(paste0("records:", t), format(cp$records[i], big.mark = "\\,"))
  put(paste0("groups:", t), format(cp$groups[i], big.mark = "\\,"))
  put(paste0("beats:", t), format(cp$beats[i], big.mark = "\\,"))
  put(paste0("p:", t), as.character(p))
  put(paste0("tail:", t), sprintf("%.1f\\%%", 100 * cp$tail10[i]))
  put(paste0("upload:", t), format(q * (p + 1), big.mark = "\\,"))
  put(paste0("covcost:", t), format(p * (p + 1) / 2, big.mark = "\\,"))
  put(paste0("factor:", t), sprintf("%.2f", p / (2 * q)))
  put(paste0("bytes:", t), format(round(q * (p + 1) * 2), big.mark = "\\,"))
}
## mean fraction of normal beats in each group of all MIT-BIH records
cl <- read.csv(out("mitbih_clusters.csv"))
nf <- sort(tapply(cl$normal, cl$cluster, mean), decreasing = TRUE)
put("normalfrac:mitbih", paste(sprintf("$%.2f$", nf), collapse = ", "))
## per-record cost of three shared bases at D = 10% in MIT-BIH
r10 <- rates[rates$corpus == "mitbih" & rates$target == 10, ]
put("extra:mitbih:10", sprintf("$%.0f$",
    r10$per_record[r10$method == "multi3"] - r10$per_record[r10$method == "spec"]))
## size of a raw MIT-BIH beat window at the 11-bit resolution of the database
put("rawbits:mitbih", as.character(11L * cp$p[cp$tag == "mitbih"]))
put("ratio:mitbih:10", sprintf("%.0f", 11 * cp$p[cp$tag == "mitbih"] /
    r10$broadcast[r10$method == "multi3"]))
## share of the distance between moment pooling and local PCA closed by three shared bases
mk <- function(t, m, k) means$beat[means$corpus == t & means$method == m & means$k == k]
for (t in tags) for (k in c(3L, 10L))
  put(sprintf("gap:%s:multi3:%d", t, k),
      pct((mk(t, "spec", k) - mk(t, "multi3", k)) / (mk(t, "spec", k) - mk(t, "local", k))))
writeLines(sprintf("\\bpval{%s}{%s}", names(val), val), out("values.tex"))

## LaTeX tables
lab <- c(local = "Local PCA (oracle)", multi3 = "Three shared bases", enorm = "Energy-normalized pooling",
         spec = "Moment pooling", covn = "Global PCA", hier = "Hierarchical alignment",
         proc = "Procrustes mean", dct = "DCT-II")
mean_of <- function(t, m, k, v = "beat") means[means$corpus == t & means$method == m & means$k == k, v]
row <- function(l, x) paste0(l, " & ", paste(x, collapse = " & "), " \\\\")
hdr <- c("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lcccccccc@{}}", "\\toprule",
         " & \\multicolumn{4}{c}{$k=3$} & \\multicolumn{4}{c}{$k=10$} \\\\",
         "\\cmidrule(lr){2-5}\\cmidrule(l){6-9}",
         " & MIT-BIH & SVDB & INCART & Finance & MIT-BIH & SVDB & INCART & Finance \\\\",
         "\\midrule")
body <- vapply(names(lab), function(m) row(lab[[m]], sprintf("%.3f",
  c(vapply(tags, mean_of, 0, m = m, k = 3), vapply(tags, mean_of, 0, m = m, k = 10)))), "")
writeLines(c(hdr, body, "\\bottomrule", "\\end{tabular*}"), out("tab_errors.tex"))

rl <- function(a, b, k, v = "beat") vapply(tags, function(t)
  ci(rules[rules$corpus == t & rules$a == a & rules$b == b & rules$k == k & rules$metric == v, ]), "")
diff_rows <- list(
  list("Three shared bases $-$ moment pooling", "multi3", "spec", 3, "beat"),
  list("\\quad ratio of traces", "multi3", "spec", 3, "trace"),
  list("\\quad $k=10$", "multi3", "spec", 10, "beat"),
  list("\\quad basis chosen from 50 beats", "multi3_n0", "spec", 3, "beat"),
  list("Energy-normalized $-$ moment pooling", "enorm", "spec", 3, "beat"),
  list("Moment pooling $-$ global PCA", "spec", "covn", 3, "beat"),
  list("Moment pooling $-$ Procrustes mean", "spec", "proc", 3, "beat"),
  list("Moment pooling $-$ hierarchical", "spec", "hier", 3, "beat"),
  list("Three shared bases $-$ hierarchical", "multi3", "hier", 3, "beat"))
body <- vapply(diff_rows, function(z) row(z[[1]], rl(z[[2]], z[[3]], z[[4]], z[[5]])), "")
writeLines(c("\\begin{tabular}{lcccc}", "\\toprule", " & MIT-BIH & SVDB & INCART & Finance \\\\",
             "\\midrule", body, "\\bottomrule", "\\end{tabular}"), out("tab_differences.tex"))

rl2 <- c(multi3 = "Three shared bases", enorm = "Energy-normalized pooling", spec = "Moment pooling",
         dct = "DCT-II", local = "Local PCA (oracle)")
fb <- function(x) ifelse(is.na(x), "--", sprintf("%.1f", x))
body <- character()
for (t in c("mitbih", "svdb", "incart")) {
  body <- c(body, sprintf("\\multicolumn{6}{@{}l}{\\emph{%s}} \\\\", CORPORA[[t]]))
  for (m in names(rl2)) {
    d <- rates[rates$corpus == t & rates$method == m, ]
    body <- c(body, paste0("\\quad ", rl2[[m]], " & ", sprintf("%.1f", d$floor[1]), " & ",
                           paste(fb(d$broadcast[match(c(20, 10, 7, 5), d$target)]), collapse = " & "), " \\\\"))
  }
}
writeLines(c("\\begin{tabular*}{\\columnwidth}{@{\\extracolsep{\\fill}}lccccc@{}}", "\\toprule",
             " & $D_{\\min}$ & $D=20$ & $10$ & $7$ & $5$ \\\\", "\\midrule", body,
             "\\bottomrule", "\\end{tabular*}"), out("tab_rates.tex"))

dl <- c(multi3 = "Three shared bases", enorm = "Energy-normalized pooling", proc = "Procrustes mean",
        local = "Local PCA (oracle)")
body <- character()
for (v in c("lda_bal", "bal")) for (k in c(3L, 10L)) {
  body <- c(body, sprintf("\\multicolumn{4}{@{}l}{\\emph{%s, $k=%d$}} \\\\",
                          if (v == "bal") "$k$-means" else "Linear discriminant", k))
  for (m in names(dl)) body <- c(body, row(paste0("\\quad ", dl[[m]]), vapply(c("mitbih", "svdb", "incart"),
    function(t) ci(down[down$corpus == t & down$metric == v & down$k == k & down$method == m, ]), "")))
}
nrec <- vapply(c("mitbih", "svdb", "incart"), function(t) down$n[down$corpus == t][1], 0)
writeLines(c("\\begin{tabular}{lccc}", "\\toprule",
             sprintf(" & MIT-BIH ($n=%d$) & SVDB ($n=%d$) & INCART ($n=%d$) \\\\", nrec[1], nrec[2], nrec[3]),
             "\\midrule", body, "\\bottomrule", "\\end{tabular}"), out("tab_downstream.tex"))

full <- do.call(rbind, lapply(tags, function(t) {
  x <- means[means$corpus == t, ]
  x[order(x$method, x$k), ]
}))
mt <- c(local = "Local PCA (oracle)", covn = "Global PCA", spec = "Moment pooling",
        enorm = "Energy-normalized pooling", hier = "Hierarchical alignment", proc = "Procrustes mean", dct = "DCT-II")
body <- character()
for (t in tags) {
  body <- c(body, sprintf("\\multicolumn{6}{@{}l}{\\emph{%s}} \\\\", CORPORA[[t]]))
  for (m in names(mt)) body <- c(body, row(paste0("\\quad ", mt[[m]]),
    sprintf("%.3f", vapply(c(1, 2, 3, 5, 10), function(k) mean_of(t, m, k), 0))))
}
writeLines(c("\\begin{tabular}{lccccc}", "\\toprule", " & $k=1$ & $2$ & $3$ & $5$ & $10$ \\\\",
             "\\midrule", body, "\\bottomrule", "\\end{tabular}"), out("tab_errors_k.tex"))
message("values: ", length(val))
