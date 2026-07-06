library(ggplot2); theme_set(theme_bw())
library(tidyverse)

dd <- readRDS("outputs/sim_batch5.rds")
print(nrow(dd))

pivfun <- function(x) {
  res <- tidyr::pivot_longer(x, contains("_"))
  if ("rvec" %in% names(res)) {
    res <- res |> mutate(rvec_f = 
                           ordered(rvec, levels = unique(rvec), labels = signif(unique(rvec), 2)))
  }
  res
}

dd_long <- pivfun(dd) |> rename(R0 = "R0vec")

gg1 <- ggplot(dd_long, aes(R0, value, colour = rvec_f)) +
  geom_point() +
  ## GAM gives nicer plots than default loess smoothing, but fails in a few cases
  ## (because quasi-eq doesn't leave enough points etc.)
  geom_smooth(method = "gam", formula = y ~ s(x)) +
  ## geom_smooth() +
  facet_wrap(~ name, scale = "free", ncol = 2) +
  scale_colour_discrete(name="r")

gg1 

ggsave("outputs/plot_batch5.pdf", width = 10, height = 5)

gg1 + scale_y_log10()
