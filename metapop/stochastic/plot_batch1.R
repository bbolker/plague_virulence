library(ggplot2); theme_set(theme_bw())
library(tidyverse)

dd0 <- readRDS("sim_batch0.rds")
dd <- readRDS("sim_batch1.rds")
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
dd_long <- pivfun(dd)

gg1 <- ggplot(dd_long, aes(R0, value, colour = alphavec_f)) +
  geom_point() +
  ## gamma = 0.75 allows slightly wigglier fits (loess is too wiggly for me)
  geom_smooth(method = "gam", formula = y ~ s(x, k = 15), method.args = list(method = "REML", gamma = 0.75)) + 
  facet_wrap(~ name, scale = "free")

gg1 + geom_point(data = dd0_long, colour =  "black", size = 3)

ggsave("plot_batch1.pdf", width = 10, height = 5)
