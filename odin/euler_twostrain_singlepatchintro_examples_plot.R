library(plagueMetapop)
library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
rot_strips <-   theme(strip.text.y = element_text(angle = 0))
zmargin <- theme(panel.spacing = grid::unit(0, "lines"))

pars <- read.csv(here::here("metapop/odin/euler_twostrain_example_pars.csv"))
dd <- readRDS(here::here("metapop/odin/outputs/euler_twostrain_singlepatchintro_examples.rds"))

m <- attr(dd, "metadata")

names(dd) <-  m$description
ddr <- dd |>
  bind_rows(.id = "parset") |>
  tidyr::pivot_longer(patches:sd, names_to = "var")

ddr_s <- filter(ddr,
                step %% 10 == 1)

oi <- palette()[-1] ## Okabe-Ito

## ugly
ggplot(ddr_s, aes(step, value, colour = state)) +
  geom_line(aes(group = interaction(run, state)), alpha = 0.3) +
  scale_colour_manual(values = oi) +
  facet_grid(parset ~ var) +
  scale_y_log10() +
  rot_strips +
  zmargin



mk_plot <- function(focal_var, log_y = FALSE) {
  d <- filter(ddr_s, var == focal_var)
  gg0 <- ggplot(ddr_s, aes(step, value, colour = state)) +
    geom_line(aes(group = interaction(run, state)), alpha = 0.3) +
    facet_wrap(~parset) +
    scale_colour_manual(values = oi) +
    labs(title = focal_var)
  if (log_y) gg0 <- gg0 + scale_y_log10()
  gg0
}

mk_plot("patches")
