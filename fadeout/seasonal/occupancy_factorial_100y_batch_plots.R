## Heatmaps of burnout and persistence for the R0 x K x alpha batch runs.
##
## Reads data/combo_summary.csv and data/settings.csv written by
## occupancy_factorial_100y_batch.R and draws one heatmap per summary
## quantity: K (x) by R0 (y), faceted scenario ~ alpha.
##
## Usage:
##   Rscript fadeout/seasonal/occupancy_factorial_100y_batch_plots.R [--pilot]
## --pilot plots the pilot run's output instead.

library(dplyr)
library(ggplot2)
library(here)

theme_set(theme_bw())

zmargin <- theme(panel.spacing = grid::unit(0, "pt"))

pilot <- "--pilot" %in% commandArgs(trailingOnly = TRUE)

outdir <- here::here(
  "fadeout", "output",
  if (pilot) "occupancy_factorial_100y_batch_pilot"
  else "occupancy_factorial_100y_batch"
)

dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)

combo_summary <- read.csv(file.path(outdir, "data", "combo_summary.csv"))

settings <- read.csv(file.path(outdir, "data", "settings.csv"))
settings <- setNames(as.list(settings$value), settings$name)


## ------------------------------------------------------------
## Heatmaps: K (x) by R0 (y), faceted scenario ~ alpha
## ------------------------------------------------------------

## Save each figure as both PDF and 300-dpi PNG (for HTML output)
save_figure <- function(plot, name, width, height) {
  for (ext in c("pdf", "png")) {
    ggsave(file.path(outdir, "figures", paste0(name, ".", ext)),
           plot, width = width, height = height, dpi = 300)
  }
}

seasonal_amp <- max(combo_summary$seasonal_amp)

scenario_labels <- c(
  seasonal = sprintf("seasonal\n(amplitude = %g)", seasonal_amp),
  constant = "constant"
)

facet_labels <- labeller(
  scenario = scenario_labels,
  alpha = function(x) paste0("alpha = ", x, "/day")
)

plot_dat <- combo_summary |>
  mutate(scenario = factor(scenario, levels = names(scenario_labels)))

heatmap <- function(var, fill_name, title, transform = "identity") {
  ggplot(plot_dat,
         aes(x = factor(K / 1000), y = factor(R0), fill = .data[[var]])) +
    geom_raster() +
    facet_grid(scenario ~ alpha, labeller = facet_labels) +
    scale_fill_viridis_c(name = fill_name, transform = transform,
                         na.value = "grey80") +
    coord_cartesian(expand = FALSE) +
    labs(
      x = "Carrying capacity per patch (K, thousands)",
      y = expression(R[0]),
      title = title,
      subtitle = sprintf(
        "%d replicate(s) per cell, %d patches, %d years",
        settings$n_reps, settings$n_patch, settings$n_years
      )
    ) +
    zmargin
}

p_occ <- heatmap(
  "occ_median",
  "median\noccupancy",
  sprintf("Patch occupancy over the last %d years", settings$summary_years)
) +
  ## Mark cells where at least one replicate went globally extinct
  geom_point(
    data = filter(plot_dat, p_persist < 1),
    shape = 4,
    colour = "white"
  ) +
  labs(caption = "x: at least one replicate went globally extinct")

p_persist <- heatmap(
  "p_persist",
  "P(global\npersistence)",
  sprintf(
    "Proportion of runs with infection persisting %d years",
    settings$n_years
  )
)

p_burnout <- heatmap(
  "p_burnout",
  "P(burnout)",
  sprintf(
    "Burnout probability after a major outbreak (extinction within %g T0)",
    settings$threshold_multiplier
  )
)

p_hazard <- heatmap(
  "fadeout_hazard",
  "fade-outs per\npatch-year",
  "Fade-out hazard after establishment",
  transform = "log10"
)

fig_width <- 16
fig_height <- 6

save_figure(p_occ, "heatmap_occupancy", fig_width, fig_height)
save_figure(p_persist, "heatmap_persistence", fig_width, fig_height)
save_figure(p_burnout, "heatmap_burnout", fig_width, fig_height)
save_figure(p_hazard, "heatmap_fadeout_hazard", fig_width, fig_height)

cat("Output: ", outdir, "\n", sep = "")
