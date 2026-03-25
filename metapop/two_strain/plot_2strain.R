library(ggplot2); theme_set(theme_bw())
library(tidyverse)

dd <- readRDS("outputs_2strain/sim_2strain.rds")
print(nrow(dd))

plot_vars <- c(
  "extinction_rate", "mean_extinct_time",
  "infected_patches", "total_inf"
)

invade_year <- attr(dd, "params0")[["invade_year"]]

dd_long <- dd %>%
  select(R0vec, alphavec, rhovec, c0vec, rvec,
         extinction_rate_1, extinction_rate_2,
         mean_extinct_time_1, mean_extinct_time_2,
         infected_patches_1, infected_patches_2,
         total_inf_1, total_inf_2) %>%
  pivot_longer(
    cols = matches("_1$|_2$"),
    names_to = c("metric", "strain"),
    names_pattern = "(.+)_([12])$"
  ) %>%
  mutate(
    strain = ifelse(strain == "1", "Resident", "Invader"),
    metric = factor(metric, levels = plot_vars),
    alphavec_f = ordered(alphavec, levels = unique(alphavec),
                         labels = signif(unique(alphavec), 2))
  )

param_combos <- dd %>%
  distinct(rhovec, c0vec, rvec) %>%
  arrange(rhovec, c0vec, rvec)

pdf("outputs_2strain/plot_2strain_points.pdf", width = 14, height = 10)

for (i in 1:nrow(param_combos)) {
  for (a in sort(unique(dd$alphavec))) {
    
    dd_sub <- dd_long %>%
      filter(rhovec == param_combos$rhovec[i],
             c0vec == param_combos$c0vec[i],
             rvec == param_combos$rvec[i],
             alphavec == a)
    
    gg <- ggplot(dd_sub, aes(x = R0vec, y = value, colour = strain)) +
      geom_point(size = 2.5, alpha = 0.8) +
      geom_hline(
        data = dd_sub %>% filter(metric == "mean_extinct_time"),
        aes(yintercept = invade_year),
        linetype = "dashed",
        colour = "black",
        linewidth = 0.6,
        inherit.aes = FALSE
      ) +
      facet_wrap(~ metric, scales = "free_y", ncol = 2,
                 labeller = as_labeller(c(
                   mean_extinct_time = "Mean Extinction Time",
                   extinction_rate = "Extinction Rate",
                   infected_patches = "Infected Patches (quasi-eq)",
                   total_inf = "Total Infections (quasi-eq)"
                 ))) +
      scale_colour_manual(values = c("Resident" = "blue", "Invader" = "red")) +
      labs(
        subtitle = paste0(
          "alpha = ", signif(a, 2),
          " | rho = ", param_combos$rhovec[i],
          ", c0 = ", param_combos$c0vec[i],
          ", r = ", param_combos$rvec[i],
          " | Invader R0 = Resident R0 - 0.2",
          " | Invader introduction year = ", invade_year
        ),
        x = "Resident R0",
        y = "Value"
      ) +
      theme(legend.position = "bottom")
    
    print(gg)
    cat("Plotted combo", i, "alpha", signif(a, 2), "\n")
  }
}

dev.off()
cat("Saved: outputs_2strain/plot_2strain_points.pdf\n")