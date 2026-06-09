library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(patchwork)
library(ggnewscale)

plot_height <- 8; plot_width <- 10
res <- readRDS(here::here("metapop/odin/outputs/ode_trough.rds"))

make_raster <- function(dat, fill_var, fill_label = fill_var,
                        title = fill_var, small_vals = NULL) {
  gg0 <- ggplot(dat, aes(x = log10(r), y = beta, fill = .data[[fill_var]])) +
    geom_raster() +
      scale_fill_viridis_c(name = fill_label, na.value = "grey80",
                           trans = "log10") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = expression(log[10]('growth')), y = expression(R[0]), title = title)

  
  if (!is.null(small_vals)) {
    ## colour 'small vals' region in pink
    gg0 <- gg0 +
      new_scale_fill() +
      geom_raster(aes(fill = factor(.data[[small_vals]]))) +
      scale_fill_manual(values = c("#00000001", "#ff0000cc"),
                        guide = "none")
  }

  return(gg0)
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
  plots <- lapply(panels, \(p) {
    small_vals <- if (p$var %in% c("Imin", "t_Imin")) "any_small" else NULL
    make_raster(dat, fill_var = p$var,
                fill_label = "",
                title =  p$label,
                small_vals = small_vals)
  })

  pw <- wrap_plots(plots, ncol = 3) +
    plot_annotation(title = dem)

  sfun <- function(ext) {
    ggsave(paste0("ode_trough_", dem, ext),
           pw, width = plot_width, height = plot_height)
  }
  sfun(".png")
  sfun(".pdf")

}
invisible(NULL)
