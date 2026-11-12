## Step 4. Figures of the article, sized for a two-column page (7 in wide).

source(file.path(Sys.getenv("BP_ARTICLE", system.file("article", package = "beatpool")), "helpers.R"))
suppressPackageStartupMessages(library(ggplot2))

lab <- c(local = "Local PCA (oracle)", covn = "Global PCA", spec = "Moment pooling",
         enorm = "Energy-normalized pooling", multi3 = "Three shared bases",
         hier = "Hierarchical alignment", proc = "Procrustes mean", dct = "DCT-II")
col <- c(local = "#000000", covn = "#7F7F7F", spec = "#E69F00", enorm = "#D55E00",
         multi3 = "#0072B2", hier = "#009E73", proc = "#56B4E9", dct = "#CC79A7")
shp <- c(local = 16, covn = 1, spec = 17, enorm = 15, multi3 = 18, hier = 4, proc = 0, dct = 3)
names(col) <- names(shp) <- lab[names(col)]
th <- theme_bp(base_size = 8) +
  theme(strip.text = element_text(face = "bold", size = 8),
        axis.title = element_text(size = 8), axis.text = element_text(size = 7),
        legend.text = element_text(size = 7.5), legend.title = element_blank(),
        legend.key.width = unit(1.4, "lines"), legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-4, 0, 0, 0), panel.spacing.x = unit(0.9, "lines"),
        plot.margin = margin(2, 6, 2, 2))
corpus_f <- function(x) factor(CORPORA[x], levels = CORPORA)
save_fig <- function(p, name, w, h) ggsave(out(paste0(name, ".pdf")), p, width = w, height = h,
                                           device = grDevices::cairo_pdf)
means <- read.csv(out("means.csv"))
means$Corpus <- corpus_f(means$corpus)

## Error against the number of retained components.
m <- means[means$method %in% c("local", "covn", "spec", "hier", "proc", "dct"), ]
m$Operator <- factor(lab[m$method], levels = lab[c("local", "covn", "spec", "hier", "proc", "dct")])
m$dash <- m$method == "covn"
p <- ggplot(m, aes(k, beat, colour = Operator, shape = Operator)) +
  geom_line(aes(linetype = dash), linewidth = 0.45, show.legend = FALSE) +
  geom_point(size = 1.4, stroke = 0.5) +
  facet_wrap(~ Corpus, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = c(1, 2, 3, 5, 10)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_colour_manual(values = col) + scale_shape_manual(values = shp) +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "22")) +
  labs(x = "Retained components, k", y = "Held-out error per beat") +
  guides(colour = guide_legend(nrow = 1)) + th
save_fig(p, "fig_error_k", 7, 2.35)

## Error against the number of shared bases, with local PCA as reference.
m <- means[means$method %in% c("enorm", "multi2", "multi3", "multi4") & means$k %in% c(3, 10), ]
m$C <- c(enorm = 1, multi2 = 2, multi3 = 3, multi4 = 4)[m$method]
m$K <- factor(paste0("k = ", m$k), levels = c("k = 3", "k = 10"))
ref <- means[means$method == "local" & means$k %in% c(3, 10), ]
ref$K <- factor(paste0("k = ", ref$k), levels = c("k = 3", "k = 10"))
p <- ggplot(m, aes(C, beat)) +
  geom_hline(data = ref, aes(yintercept = beat, linetype = "Local PCA (oracle)"),
             colour = "grey35", linewidth = 0.4) +
  geom_line(aes(colour = "Shared bases from transmitted moments"), linewidth = 0.5) +
  geom_point(aes(colour = "Shared bases from transmitted moments"), size = 1.4) +
  facet_grid(K ~ Corpus, scales = "free_y") +
  scale_colour_manual(values = "#0072B2") +
  scale_linetype_manual(values = "22") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
  labs(x = "Number of shared bases, C (C = 1 is energy-normalized pooling)", y = "Held-out error per beat") + th
save_fig(p, "fig_shared_bases", 7, 3.1)

## One point per record.
h <- do.call(rbind, lapply(c("mitbih", "svdb", "incart"), function(t) {
  d <- read.csv(out(sprintf("%s_heldout.csv", t)), colClasses = c(record = "character"))
  d <- d[d$k == 3 & d$method %in% c("spec", "multi3"), ]
  w <- reshape(d[, c("record", "method", "beat")], idvar = "record", timevar = "method",
               direction = "wide")
  data.frame(Corpus = corpus_f(t), spec = w$beat.spec, multi3 = w$beat.multi3, row.names = NULL)
}))
lim <- c(0, ceiling(max(h$spec, h$multi3) * 10) / 10)
p <- ggplot(h, aes(spec, multi3)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey45", linewidth = 0.4) +
  geom_point(colour = "#0072B2", size = 1.1, alpha = 0.8) +
  facet_wrap(~ Corpus, nrow = 1) +
  coord_equal(xlim = lim, ylim = lim, expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) + scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  labs(x = "Moment pooling", y = "Three shared bases") + th
save_fig(p, "fig_records", 7, 2.55)

## Rate--distortion with one broadcast of the shared bases.
rc <- do.call(rbind, lapply(c("mitbih", "svdb", "incart"), function(t) {
  d <- read.csv(out(sprintf("%s_rate_curve.csv", t)))
  d$Corpus <- corpus_f(t)
  d
}))
rc <- rc[rc$method %in% c("local", "spec", "enorm", "multi3", "dct"), ]
rc$Operator <- factor(lab[rc$method], levels = lab[c("local", "multi3", "enorm", "spec", "dct")])
p <- ggplot(rc[order(rc$broadcast), ], aes(broadcast, D, colour = Operator, shape = Operator)) +
  geom_path(linewidth = 0.45) + geom_point(size = 1.1) +
  facet_wrap(~ Corpus, nrow = 1) +
  coord_cartesian(xlim = c(15, 75), ylim = c(0, 30), expand = FALSE) +
  scale_colour_manual(values = col) + scale_shape_manual(values = shp) +
  labs(x = "Bits per beat", y = "Distortion D (%)") +
  guides(colour = guide_legend(nrow = 1)) + th
save_fig(p, "fig_rates", 7, 2.45)
message("figures written to ", out_dir)
