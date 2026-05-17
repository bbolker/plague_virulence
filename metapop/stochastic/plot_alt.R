#!/usr/bin/env Rscript
library(optparse)
op.parser <- OptionParser(prog="sim_onestrain_poisson",
                          option_list = list(
                            make_option(c("-i", "--input"), "store",
                                        help = "input data file",
                                        default = "outputs/sim.rds"),
                            make_option(c("-o", "--output"), "store",
                                        help = "output plot file",
                                        default = "outputs/poisson_plot_alt.pdf"),
                            make_option(c("-s", "--smooth"), "store",
                                        help = "smoother type (options: 'none', 'loess', 'line')",
                                        default = "loess"),
                            make_option(c("-w", "--width"), "store",
                                        help = "plot width",
                                        default = 20),
                            make_option(c("--height"), "store",
                                        help = "plot height",
                                        default = 10),
                            make_option(c("--bw", "store_false"),
                                        help = "black and white plot",
                                        default = FALSE),
                            make_option(c("--which", "store",
                                          help = "plots to create (default = NA)",
                                          default = NA_integer_))
                          )
                          )
                            
opt <- parse_args(op.parser)

opt$input <- "outputs/sim_batch1.rds"

library(tidyverse); theme_set(theme_bw())
zmargin <- theme(panel.spacing = grid::unit(0, "lines"))
library(colorspace)
library(cowplot)

dd <- readRDS(opt$input) ##  |> na.omit()

plot_vars <- c("extinction_rate", "mean_extinct_time",
               "total_pops", "infected_patches", "total_inf", "total_pops_qe")

id_vars <- c("R0vec", "alphavec", "rhovec", "c0vec", "rvec")
dd_long <- dd |>
  mutate(across(c(total_pops, total_inf), ~ . / 1e6)) |>
  select(any_of(id_vars), any_of(plot_vars)) |>
  pivot_longer(cols = any_of(plot_vars),
               names_to = "metric") |>
  mutate(
    metric = factor(metric, levels = plot_vars),
    log10alpha = round(log10(alphavec),1),
    alphavec_f = factor(alphavec, levels = unique(alphavec),
                        labels = signif(unique(alphavec), 2))
  ) |>
  na.omit()

## quasi-binomial smooth (slow)
qb_smooth <- geom_smooth(method = "gam",
                         method.args = list(family = quasibinomial),
                         alpha = 0.2)

## quasi-poisson smooth (a little less slow)
qp_smooth <- geom_smooth(method = "gam",
                         method.args = list(family = quasipoisson),
                         alpha = 0.2)

## heat_hcl args
## (n, h = c(0, 90), c. = c(100, 30), l = c(50, 90), power = c(1/5, 
##     1), gamma = NULL, fixup = TRUE, alpha = 1, ...) 

hpal <-   scale_color_continuous_sequential(
  l1 = 50, l2 = 70, h1 = 0, h2 = 90, c1 = 100, c2 = 30,
  guide = guide_legend())
  # use guide_legend rather than guide_colourbar since we have discrete values anyway


## use log10-alpha (prettier for non-integer values)
gg0 <- ggplot(dd_long,
              ## aes(R0vec, value, colour = log10alpha, group = log10alpha)) +
              aes(R0vec, value, group = log10alpha)) +
  geom_point() +
  hpal 

if (!opt$bw) gg0 <- gg0 + aes(colour = log10alpha)
## if we did want a reversed colour bar (to match facet ordering top-to-bottom:
##  guide_colorbar(reverse=TRUE)
## https://aosmith.rbind.io/2018/01/19/reversing-the-order-of-a-ggplot2-legend/


gg0 + geom_line() + facet_grid(metric ~ log10alpha, scale = "free",
                               labeller = label_both) +
  zmargin +
    theme(strip.text.y = element_text(angle = 0)) +
    guides(color = "none") +
    labs(x = expression(R[0]))

ggsave("burnout_sum.png", width = 10, height = 8)
