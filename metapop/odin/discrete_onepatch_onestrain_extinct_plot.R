library(dplyr)
library(ggplot2)
library(patchwork)

theme_set(theme_bw())

dat <- readRDS("discrete_onepatch_onestrain_extinct.rds") |>
  mutate(log10K = log10(K))

make_raster <- function(fill_var, fill_label) {
  ggplot(dat, aes(R0, log10K, fill = .data[[fill_var]])) +
    geom_raster() +
    scale_fill_viridis_c(trans = "log10", na.value = "grey80",
                         name = fill_label) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = 3:6,
                       labels = scales::label_math(10^.x)) +
    labs(x = expression(R[0]), y = expression(log[10](K)))
}

p1 <- make_raster("ext_prob.I1",      "extinction\nprobability")
p2 <- make_raster("mean_ext_time.I1", "mean extinction\ntime (steps)")

print(p1 + p2)
ggsave("discrete_onepatch_onestrain_extinct.png", width = 10, height = 5)
