library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(patchwork)

plot_height <- 8; plot_width <- 10
res <- readRDS(here::here("odin/outputs/ode_trough.rds"))

res2 <- res |>
  mutate(across(c(t_Imin, Imin),
                ~ifelse(!is.finite(any_small) | any_small, NaN, .)))

make_raster <- function(dat, fill_var, fill_label = fill_var,
                        title = fill_var, small_vals = NULL,
                        contour = TRUE, zlimits = c(NA, NA)) {
  gg0 <- ggplot(dat, aes(x = log10(r), y = beta)) +
    geom_raster(aes(fill = .data[[fill_var]])) +
      scale_fill_viridis_c(name = fill_label, na.value = "grey80",
                           trans = "log10",
                           limits = zlimits) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = expression(log[10]('growth')), y = expression(R[0]), title = title)

  if (contour) {
      gg0 <- gg0 + geom_contour(aes(z = log10(.data[[fill_var]])), colour = "red")
  }

  if (length(unique(dat$demography)) > 1)
    gg0 <- gg0 + facet_wrap(~demography, nrow = 1)

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
  dat <- filter(res2, demography == dem)
  plots <- lapply(panels, \(p) {
    small_vals <- if (p$var %in% c("Imin", "t_Imin")) "any_small" else NULL
    make_raster(dat, fill_var = p$var,
                fill_label = "",
                title =  p$label,
                small_vals = small_vals)
  })

  pw <- wrap_plots(plots, ncol = 3) +
    plot_annotation(title = dem)

  sfun <- function(ext, base = "ode_trough_", plot = pw) {
    ggsave(paste0(base, dem, ext),
           plot, width = plot_width, height = plot_height)
  }
  sfun(".png")
  sfun(".pdf")

}

## Combined plot: both demographies, Imin and trough_area only
## Top row: two separate Imin plots with independent scale bars
Imin_plots <- lapply(c("logistic", "linear"), \(dem) {
  res_filtered <- dplyr::filter(res2, demography == dem)
  make_raster(res_filtered,
              fill_var    = "Imin",
              fill_label  = "",
              title       = paste("I_min", dem),
              small_vals  = "any_small",
              zlimits     = c(max(1e-100, min(res_filtered$Imin, na.rm = TRUE)),
                                  NA))
})

## Bottom row: trough_area faceted by demography (unchanged)
trough_plot <- make_raster(res, fill_var   = "trough_area",
                               fill_label = "",
                               title      = "trough area",
                               zlimits    = c(NA, NA))

pw_combined <- (Imin_plots[[1]] | Imin_plots[[2]]) / trough_plot

plot_height <- 8; plot_width <- 10
dem <- ""
sfun(".png", base = "ode_trough_combined", plot = pw_combined)
sfun(".pdf", base = "ode_trough_combined", plot = pw_combined)

ggplot(res, aes(Imin, trough_area)) +
    geom_point(aes(colour = demography, shape = demography)) +
    scale_colour_brewer(palette = "Dark2") +
    scale_x_log10() + scale_y_log10() +
    facet_wrap(~demography, scale = "free")
ggsave("Imin_vs_trough.pdf", width = 12, height = 8)
