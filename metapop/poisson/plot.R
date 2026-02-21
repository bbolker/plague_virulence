library(ggplot2); theme_set(theme_bw())
library(tidyverse)

dd <- readRDS("outputs/sim.rds")
print(nrow(dd))

plot_vars <- c("extinction_rate", "mean_extinct_time",
               "total_pops", "infected_patches", "total_inf")

dd_long <- dd %>%
  select(R0vec, kappavec, etavec, c0vec, rvec, all_of(plot_vars)) %>%
  pivot_longer(cols = all_of(plot_vars), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = plot_vars),  
    kappavec_f = ordered(kappavec, levels = unique(kappavec),
                         labels = signif(unique(kappavec), 2))
  )

param_combos <- dd %>%
  distinct(etavec, c0vec, rvec) %>%
  arrange(etavec, c0vec, rvec)

pdf("outputs/plot.pdf", width = 14, height = 10)

for (i in 1:nrow(param_combos)) {
  
  dd_sub <- dd_long %>%
    filter(etavec == param_combos$etavec[i],
           c0vec == param_combos$c0vec[i],
           rvec == param_combos$rvec[i])
  
  gg <- ggplot(dd_sub, aes(R0vec, value, colour = kappavec_f)) +
    geom_line(linewidth = 1, alpha = 0.7) +
    geom_point(size = 2.5) +
    facet_wrap(~ metric, scales = "free_y", ncol = 3,
               labeller = as_labeller(c(
                 mean_extinct_time = "Mean Extinction Time",
                 extinction_rate = "Extinction Rate",
                 total_pops = "Total Population (quasi-eq)",
                 infected_patches = "Infected Patches (quasi-eq)",
                 total_inf = "Total Infections (quasi-eq)"
               ))) +
    scale_colour_discrete(name = "kappa") +
    labs(subtitle = paste0("eta = ", param_combos$etavec[i], 
                           ",  c0 = ", param_combos$c0vec[i],
                           ",  r = ", param_combos$rvec[i]),
         x = "R0", y = "Value") +
    theme(legend.position = "bottom")
  
  print(gg)
  
  cat("Plotted combo", i, "of", nrow(param_combos), "\n")
}

dev.off()

cat("Saved: outputs/plot.pdf\n")