library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(patchwork)

res <- readRDS("outputs/ode_trough.rds")

make_raster <- function(dat, fill_var, fill_label = fill_var, title = fill_var) {
  ggplot(dat, aes(x = log10(r), y = beta, fill = .data[[fill_var]])) +
    geom_raster() +
      scale_fill_viridis_c(name = fill_label, na.value = "grey80",
                           trans = "log10") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = expression(log[10]('growth')), y = expression(R[0]), title = title)
}

panels <- list(
  list(var = "eq_S",             label = "equil S*"),
  list(var = "eq_I",             label = "equil I*"),
  list(var = "t_enter.boundary", label = "t(enter boundary)"),
  list(var = "t_Imin",           label = "t(I_min)"),
  list(var = "Imin",             label = "I_min"),
  list(var = "t_Smin",           label = "t(S_min)"),
  list(var = "Smin",             label = "S_min"),
  list(var = "t_leave.boundary", label = "t(leave boundary)"),
  list(var = "trough_area",      label = "trough area")
)

for (dem in c("logistic", "linear")) {
  dat <- filter(res, demography == dem)
  plots <- lapply(panels, \(p) make_raster(dat, fill_var = p$var,
                                           fill_label = "",
                                           title =  p$label))
  pw <- wrap_plots(plots, ncol = 3) +
    plot_annotation(title = dem)
  ggsave(paste0("ode_trough_", dem, ".png"), pw, width = 16, height = 8)
  ggsave(paste0("ode_trough_", dem, ".pdf"), pw, width = 16, height = 8)
}
invisible(NULL)
