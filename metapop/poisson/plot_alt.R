#!/usr/bin/env Rscript
library(optparse)
op.parser <- OptionParser(prog="sim_onestrain_poisson",
                          option_list = list(
                            make_option(c("-i", "--input"), "store",
                                        help = "input data file",
                                        default = "outputs/sim.rds"),
                            make_option(c("-o", "--output"), "store",
                                        help = "output plot file",
                                        default = "outputs/poisson_plot_alt.pdf"))
                          )
                            
opt <- parse_args(op.parser)

library(tidyverse); theme_set(theme_bw())
zmargin <- theme(panel.spacing = grid::unit(0, "lines"))
library(colorspace)
library(cowplot)

dd <- readRDS(opt$input)

plot_vars <- c("extinction_rate", "mean_extinct_time",
               "total_pops", "infected_patches", "total_inf")

dd_long <- dd |>
  mutate(across(c(total_pops, total_inf), ~ . / 1e6)) |>
  select(R0vec, alphavec, rhovec, c0vec, rvec, all_of(plot_vars)) |>
  pivot_longer(cols = all_of(plot_vars),
               names_to = "metric") |>
  mutate(
    metric = factor(metric, levels = plot_vars),
    log10alpha = round(log10(alphavec),1),
    alphavec_f = factor(alphavec, levels = unique(alphavec),
                        labels = signif(unique(alphavec), 2))
  )

## quasi-binomial smooth (slow)
qb_smooth <- geom_smooth(method = "gam",
                         method.args = list(family = quasibinomial),
                         alpha = 0.2)

## quasi-poisson smooth (a little less slow)
qp_smooth <- geom_smooth(method = "gam",
                         method.args = list(family = quasipoisson),
                         alpha = 0.2)

## use log10-alpha (prettier for non-integer values)
gg0 <- ggplot(dd_long,
              ## aes(R0vec, value, colour = log10alpha, group = log10alpha)) +
              aes(R0vec, value, group = log10alpha)) +
  geom_point() +
  # use guide_legend rather than guide_colourbar since we have discrete values anyway
  scale_color_continuous_sequential(palette = "Heat",
                                    guide = guide_legend())

## if we did want a reversed colour bar (to match facet ordering top-to-bottom:
##  guide_colorbar(reverse=TRUE)
## https://aosmith.rbind.io/2018/01/19/reversing-the-order-of-a-ggplot2-legend/

plot_fun <- function(focal_metric = "extinction_rate", add_smooth = FALSE, limits = c(0, NA),
                     title = focal_metric) {
  dd2 <- dd_long |> filter(metric  == focal_metric,
                           ## restrict to 'interesting' cases
                           log10alpha>=(-5.5),
                           log10alpha<=(-3.5))

  rho_labs <- function(values) {
    lapply(labels, function(values) {
      browser()
      values <- paste0("list(",
                       sprintf("rho: %s", values),
                       ")")
      lapply(values, function(expr) c(parse(text = expr)))
    })
  }
  
  gg1 <- gg0 + dd2 +
    geom_line() +
    facet_grid(log10alpha~rhovec,
               ## https://stackoverflow.com/a/74698645/190277`
               labeller = labeller(
                 rhovec = as_labeller(~paste0("rho: ", .x), label_parsed),
                 log10alpha = as_labeller(~paste0("log[10](alpha): ", .x), label_parsed))) +
    scale_y_continuous(limits = limits,
                       oob = scales::squish) +
    labs(title = title) +
    zmargin

  if (add_smooth) {
    ## https://stackoverflow.com/a/46285325/190277
    gg1 <- gg1 + stat_smooth(geom = "line", lty = 2, alpha = 0.5, se = FALSE)
  }
  return(gg1)
  
}

design <- tribble(
  ~ focal_metric, ~add_smooth, ~title,
  "mean_extinct_time", FALSE, "Mean Extinction Time",
  "extinction_rate", FALSE, "Extinction Rate",
  "total_pops", TRUE,  "Total host population (x 1e6; quasi-eq)",
  "infected_patches", TRUE, "Infected patches (quasi-eq)",
  "total_inf", TRUE, "Total Infections (x 1e6; quasi-eq)")

plot_list <- list()
for (i in 1:nrow(design)) {
  plot_list[[design$focal_metric[i]]] <-
    do.call(plot_fun, design[i,])
}

if (FALSE) do.call(plot_fun, design[1,])

## cleanup: colour legend only on last plot
for (i in 1:(length(plot_list)-1)) {
  plot_list[[i]] <- plot_list[[i]] + guides(color = "none")
}

## cleanup: no y-axis labels (should apply upstream ...)
plot_list <- lapply(plot_list, \(x) (x + labs(y="")))

pdf(opt$output, width = 20, height = 10)
suppressWarnings(
  plot_grid(plotlist = plot_list)
)
dev.off()
