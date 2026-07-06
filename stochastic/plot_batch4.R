library(ggplot2); theme_set(theme_bw())
library(tidyverse)

dd <- readRDS("outputs/sim_batch4.rds")
print(nrow(dd))

pivfun <- function(x) {
  res <- tidyr::pivot_longer(x, contains("_"))
  if ("c0vec" %in% names(res)) {
    res <- res |> mutate(c0vec_f = 
                           ordered(c0vec, levels = unique(c0vec), labels = signif(unique(c0vec), 2)))
  }
  res
}

dd_long <- pivfun(dd) |> rename(R0 = "R0vec")

gg1 <- ggplot(dd_long, aes(R0, value, colour = c0vec_f)) +
  geom_point() +
  ## GAM gives nicer plots than default loess smoothing, but fails in a few cases
  ## (because quasi-eq doesn't leave enough points etc.)
  geom_smooth(method = "gam", formula = y ~ s(x)) +
  ## geom_smooth() +
  facet_wrap(~ name, scale = "free", ncol = 2) +
  scale_colour_discrete(name="c0")

gg1 

ggsave("outputs/plot_batch4.pdf", width = 10, height = 5)

gg1 + scale_y_log10()
