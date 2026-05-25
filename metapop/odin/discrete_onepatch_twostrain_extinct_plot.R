library(ggplot2); theme_set(theme_bw())

## temporary, from checkpoint


fn <- "discrete_onepatch_twostrain_extinct_checkpoint.rds"
R0vec <- seq(1, 5, by = 0.025)
dd    <- expand.grid(R01 = R0vec, R02 = R0vec)
res <- readRDS(fn)
res_df <- t(res) |> as.data.frame() |> dplyr::bind_cols(dd)
dd2    <- res_df |> tidyr::pivot_longer(cols = c("I1", "I2"))
dd2_i1 <- dd2 |> dplyr::filter(name == "I1")


## dd2_i1 <- readRDS(fn)

ggplot(dd2_i1, aes(R01, R02, fill = value)) +
  geom_raster() +
    scale_fill_viridis_c(trans = "log10", na.value = 
                         adjustcolor("yellow", alpha.f = 0.5)
                         ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  labs(x = expression(R[0]*' of strain 1'),
       y = expression(R[0]*' of strain 2'),
       title = expression('extinction time of strain 1 '*(list(K==10^6, r==0.125)))) +
  annotate(x = 2.25, y = 1.5, size = 5, label = "time > 1000", geom = "label")

ggsave(width = 6.5, height = 6, "discrete_onepatch_twostrain_extinct.png")
