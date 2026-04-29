library(ggplot2); theme_set(theme_bw())
library(tidyverse)
library(colorspace)

dd <- readRDS("outputs/sim.rds")
print(nrow(dd))

plot_vars <- c("extinction_rate", "mean_extinct_time",
               "total_pops", "infected_patches", "total_inf")

dd_long <- dd %>%
  select(R0vec, alphavec, rhovec, c0vec, rvec, all_of(plot_vars)) %>%
  pivot_longer(cols = all_of(plot_vars), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = plot_vars),  
    alphavec_f = factor(alphavec, levels = unique(alphavec),
                        labels = signif(unique(alphavec), 2))
  )

param_combos <- dd %>%
  distinct(rhovec, c0vec, rvec) %>%
  arrange(rhovec, c0vec, rvec)

pdf("outputs/plot.pdf", width = 14, height = 10)

qb_smooth <- geom_smooth(method = "gam",
                         method.args = list(family = quasibinomial),
                         alpha = 0.2)

qp_smooth <- geom_smooth(method = "gam",
                         method.args = list(family = quasipoisson),
                         alpha = 0.2)

gg0 <- ggplot(dd_long,
              aes(R0vec, value, colour = alphavec, group = alphavec)) +
  geom_point() +
  scale_color_continuous_sequential(palette = "Heat",
                                    trans = "log10",
                                    breaks = unique(dd_long$alphavec),
                                    labels = round(log10(unique(
                                      dd_long$alphavec)), 1))
               
gg0 + (dd_long |> filter(metric  == "extinction_rate")) +
  ##  qb_smooth +
  ## geom_smooth(method = "loess", method.args = list(span  = 0.01)) +
  geom_line() +
  facet_grid(alphavec~rhovec, labeller = label_both) +
  scale_y_continuous(limits = c(0,1),
                     oob = scales::squish) +
  labs(title = "extinction rate")

gg0 + (dd_long |> filter(metric  == "infected_patches")) +
  geom_smooth() +
  facet_grid(alphavec~rhovec, labeller = label_both) +
  scale_y_continuous(limits=c(0,NA), oob = scales::squish)


