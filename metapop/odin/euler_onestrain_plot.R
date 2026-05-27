#!/usr/bin/env Rscript
library(optparse)

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-i", "--input"), default = "sharcnet/outputs/euler_onestrain.rds",
              help = "input RDS file [default %default]"),
  make_option(c("-o", "--output"), default = "euler_onestrain_plot.pdf",
              help = "output plot file [default %default]"),
  make_option(c("-w", "--width"),  type = "double", default = 14,
              help = "plot width in inches [default %default]"),
  make_option(c("--height"), type = "double", default = 10,
              help = "plot height in inches [default %default]")
)))

library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(tidyr)
library(cowplot)
library(colorspace)

zmargin <- theme(panel.spacing = grid::unit(0, "lines"))

batch2_file <- sub("\\.rds$", "_batch2.rds", opt$input)
dd <- dplyr::bind_rows(
  readRDS(opt$input),
  if (file.exists(batch2_file)) readRDS(batch2_file) else NULL
)

plot_vars <- c("ext_prob.I1", "mean_ext_time.I1",
               "qe_infpop.I1", "qe_infpatch.I1", "qe_pop_S.I1")

dd_long <- dd |>
  mutate(log10K     = round(log10(K), 2),
         log10alpha = round(log10(alpha), 1)) |>
  select(R0, log10K, log10alpha, any_of(plot_vars)) |>
  pivot_longer(cols = any_of(plot_vars), names_to = "metric") |>
  mutate(metric   = factor(metric, levels = plot_vars),
         log10K_f = factor(log10K))

hpal <- scale_color_continuous_sequential(
  l1 = 50, l2 = 70, h1 = 0, h2 = 60, c1 = 80, c2 = 30,
  name = expression(log[10](K)),
  guide = guide_legend())

my_facet <- facet_wrap(
  ~ log10alpha,
  labeller = labeller(
    log10alpha = as_labeller(~ paste0("log[10](alpha): ", .x), label_parsed)))

plot_fun <- function(focal_metric, title = focal_metric) {
  dd2 <- filter(dd_long, metric == focal_metric, !is.na(value))
  ggplot(dd2, aes(R0, value, colour = log10K, group = log10K_f)) +
    geom_line() +
    geom_point(size = 0.8) +
    my_facet +
    hpal +
    labs(title = title, x = expression(R[0]), y = "") +
    zmargin
}

design <- tribble(
  ~focal_metric,      ~title,
  "ext_prob.I1",      "Extinction probability",
  "mean_ext_time.I1", "Mean extinction time (steps)",
  "qe_infpop.I1",     "QE infected population",
  "qe_infpatch.I1",   "QE infected patches",
  "qe_pop_S.I1",      "QE susceptible population"
)

has_data <- dd_long |> filter(!is.na(value)) |> pull(metric) |> unique() |> as.character()
design <- filter(design, focal_metric %in% has_data)
plot_list <- Map(plot_fun, design$focal_metric, design$title)

## colour legend only on last panel
for (i in seq_len(length(plot_list) - 1)) {
  plot_list[[i]] <- plot_list[[i]] + guides(color = "none")
}

g <- suppressWarnings(plot_grid(plotlist = plot_list))
for (ext in c("pdf", "png")) {
  ggsave(sub("\\.[^.]+$", paste0(".", ext), opt$output), g,
         width = opt$width, height = opt$height)
}

