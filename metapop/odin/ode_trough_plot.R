library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(patchwork)

res <- readRDS("outputs/ode_trough.rds")

make_raster <- function(dat, fill_var, fill_label = fill_var, title = fill_var) {
  ggplot(dat, aes(x = log10(r), y = beta, fill = .data[[fill_var]])) +
    geom_raster() +
    scale_fill_viridis_c(name = fill_label, na.value = "grey80") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = expression(log[10](r)), y = expression(beta), title = title)
}

panels <- list(
  list(var = "eq",               label = "equilibrium\nI*"),
  list(var = "t_enter.boundary", label = "t enter\nboundary"),
  list(var = "t_Imin",           label = "t I min"),
  list(var = "Imin",             label = "I min"),
  list(var = "t_Smin",           label = "t S min"),
  list(var = "Smin",             label = "S min"),
  list(var = "t_leave.boundary", label = "t leave\nboundary"),
  list(var = "trough_area",      label = "trough\narea")
)

for (dem in c("logistic", "linear")) {
  dat <- filter(res, demography == dem)
  plots <- lapply(panels, \(p) make_raster(dat, p$var, p$label))
  pw <- wrap_plots(plots, ncol = 4) +
    plot_annotation(title = dem)
  ggsave(paste0("ode_trough_", dem, ".png"), pw, width = 16, height = 8)
  ggsave(paste0("ode_trough_", dem, ".pdf"), pw, width = 16, height = 8)
}
