library(ggplot2); theme_set(theme_bw())
library(tidyverse)

dd <- expand.grid(R0 = seq(1.5, 3, by = 0.1),
                  alphavec = 5*10^seq(-6,-3, by = 0.5))
res <- readRDS("sim_batch1.rds") |> do.call(what = "rbind")
dd <- data.frame(dd[1:nrow(res),], res)
print(nrow(res))

dd_long <- tidyr::pivot_longer(dd, -(1:2)) |>
  mutate(alphavec_f = 
           ordered(alphavec, levels = unique(alphavec), labels = signif(unique(alphavec), 2)))

ggplot(dd_long, aes(R0, value, colour = alphavec_f)) +
  geom_point() +
  geom_smooth() + 
  facet_wrap(~ name, scale = "free")

