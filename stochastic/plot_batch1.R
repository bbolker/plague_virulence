library(ggplot2); theme_set(theme_bw())
library(tidyverse)

dd0 <- readRDS("outputs/sim_batch0.rds")
dd <- readRDS("outputs/sim_batch1.rds")
print(nrow(dd))

pivfun <- function(x) {
  res <- tidyr::pivot_longer(x, contains("_"))
  if ("alphavec" %in% names(res)) {
    res <- res |> mutate(alphavec_f = 
                           ordered(alphavec, levels = unique(alphavec), labels = signif(unique(alphavec), 2)))
  }
  res
}

dd0_long <- pivfun(dd0)
dd_long <- pivfun(dd) |> rename(R0 = "R0vec")

gg1 <- ggplot(dd_long, aes(R0, value, colour = alphavec_f)) +
  geom_point() +
  ## GAM gives nicer plots than default loess smoothing, but fails in a few cases
  ## (because quasi-eq doesn't leave enough points etc.)
  geom_smooth(method = "gam", formula = y ~ s(x)) +
  ## geom_smooth() +
  facet_wrap(~ name, scale = "free", ncol = 2) +
  scale_colour_discrete(name="alpha")

gg1 + geom_point(data = dd0_long, colour =  "black", size = 3)

ggsave("plot_batch1.pdf", width = 10, height = 5)

gg1 + scale_y_log10()
