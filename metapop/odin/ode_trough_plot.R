library(ggplot2); theme_set(theme_bw())
library(dplyr)

res <- readRDS("ode_trough.rds")
res_long <- tidyr::pivot_longer(res, -c(demography:r), names_to = "var")

## want free scales on rasters, will have to patchwork etc.
res_long |> filter(demography == "linear") |>
    ggplot(aes(x = log10(r), y = beta)) +
    geom_raster(aes(fill = value)) +
    facet_wrap(~var, scale = "free") +
    scale_fill_viridis_c()
