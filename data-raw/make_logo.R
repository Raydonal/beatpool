# Hex sticker for beatpool.
# Run from the package root: Rscript data-raw/make_logo.R
library(ggplot2)
library(hexSticker)

beat_file <- c("../analysis/results/beats_mitdb.rds",
               "analysis/results/beats_mitdb.rds")
beat_file <- beat_file[file.exists(beat_file)][1]

if (!is.na(beat_file)) {
  X <- readRDS(beat_file)[[1]]
  X <- X[seq_len(min(6, nrow(X))), , drop = FALSE]
  time <- seq_len(ncol(X)) - (ncol(X) + 1) / 2
  df <- data.frame(
    time = rep(time, times = nrow(X)),
    value = as.vector(t(X)),
    beat = factor(rep(seq_len(nrow(X)), each = ncol(X)))
  )
} else {
  time <- seq(-30, 30)
  template <- -0.15 * exp(-(time / 12)^2) +
    1.1 * exp(-(time / 2)^2) - 0.45 * exp(-((time - 5) / 3)^2) +
    0.25 * exp(-((time - 18) / 7)^2)
  df <- do.call(rbind, lapply(seq_len(6), function(b) {
    data.frame(time = time, value = template * (1 + 0.08 * b) +
                 0.03 * sin(time / 2 + b), beat = factor(b))
  }))
}

img <- ggplot(df, aes(time, value, group = beat)) +
  geom_line(colour = "#8FD3C7", linewidth = 0.55, alpha = 0.55) +
  geom_line(data = df[df$beat == "1", ], colour = "#FFFFFF",
            linewidth = 0.9) +
  coord_cartesian(expand = FALSE) +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA))

dir.create("inst/figures", recursive = TRUE, showWarnings = FALSE)

hexSticker::sticker(
  img,
  package = "beatpool",
  p_size = 25, p_x = 1, p_y = 1.44, p_color = "#FFFFFF",
  p_family = "sans", p_fontface = "bold",
  s_x = 1, s_y = 0.82, s_width = 1.22, s_height = 0.82,
  h_fill = "#0B4F6C", h_color = "#F2A541",
  filename = "inst/figures/beatpool.png", dpi = 600
)

dir.create("man/figures", showWarnings = FALSE)
file.copy("inst/figures/beatpool.png", "man/figures/logo.png",
          overwrite = TRUE)
