library(ggplot2); theme_set(theme_bw())
library(tidyverse)

dd <- readRDS("outputs/sim.rds")
print(nrow(dd))

pivfun <- function(x) {
  res <- tidyr::pivot_longer(x, contains("_"))
  if ("kappavec" %in% names(res)) {
    res <- res |> mutate(kappavec_f = 
                           ordered(kappavec, levels = unique(kappavec), 
                                   labels = signif(unique(kappavec), 2)))
  }
  res
}

dd_long <- pivfun(dd) |> rename(R0 = "R0vec")

param_combos <- dd_long |> 
  distinct(etavec, c0vec, rvec) |> 
  arrange(etavec, c0vec, rvec)

pdf("outputs/plot.pdf", width = 10, height = 5)

for (i in 1:nrow(param_combos)) {
  
  dd_subset <- dd_long |> 
    filter(etavec == param_combos$etavec[i],
           c0vec == param_combos$c0vec[i],
           rvec == param_combos$rvec[i])
  
  gg1 <- ggplot(dd_subset, aes(R0, value, colour = kappavec_f)) +
    geom_point() +
    geom_smooth(method = "gam", formula = y ~ s(x)) +
    facet_wrap(~ name, scale = "free", ncol = 2) +
    scale_colour_discrete(name="kappa") +
    labs(subtitle = paste0("eta = ", param_combos$etavec[i], 
                           ",  c0 = ", param_combos$c0vec[i],
                           ",  r = ", param_combos$rvec[i]))
  
  print(gg1)
}

dev.off()

cat("Saved: outputs/plot.pdf\n")