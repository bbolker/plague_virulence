## Heatmaps of burnout and persistence for the R0 x K x alpha batch runs.
##
## Reads data/combo_summary.csv, data/run_summary.csv and data/settings.csv
## written by occupancy_factorial_100y_batch.R and draws one heatmap per
## summary quantity: K (x) by R0 (y), faceted seasonal amplitude ~ alpha, plus
## amplitude (x) by R0 (y) heatmaps faceted by K at a single alpha.
##
## Colours run from white to red. Probabilities are on a logit scale, the
## fade-out hazard on a log scale. Zeros (and ones) cannot be shown on these
## scales, so each quantity is clamped to [lo, 1 - lo] (hazard: >= lo), where
## lo is the value one event gives in the cell with the largest denominator
## (pooled over replicates), or half the smallest observed nonzero value if
## that is smaller.
##
## Usage:
##   Rscript fadeout/seasonal/occupancy_factorial_100y_batch_plots.R [--pilot]
## --pilot plots the pilot run's output instead.

library(dplyr)
library(ggplot2)
library(here)

theme_set(theme_bw())

## Facet columns abut, facet rows are separated
zmargin <- theme(
  panel.spacing.x = grid::unit(0, "pt"),
  panel.spacing.y = grid::unit(6, "pt")
)

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
## Heatmaps: K (x) by R0 (y), faceted seasonal amplitude ~ alpha
## ------------------------------------------------------------

## Save each figure as both PDF and 300-dpi PNG (for HTML output)
save_figure <- function(plot, name, width, height) {
  for (ext in c("pdf", "png")) {
    ggsave(file.path(outdir, "figures", paste0(name, ".", ext)),
           plot, width = width, height = height, dpi = 300)
  }
}

amp_label <- function(x) {
  ifelse(as.numeric(x) == 0, "constant", paste0("amplitude = ", x))
}

facet_labels <- labeller(
  seasonal_amp = amp_label,
  alpha = function(x) paste0("alpha = ", x, "/day")
)

## Pooled denominators per cell, from the per-run counts
denoms <- read.csv(file.path(outdir, "data", "run_summary.csv")) |>
  group_by(R0, K, alpha, seasonal_amp) |>
  summarise(
    n_runs = n(),
    n_post_major = sum(n_burnout + n_survived),
    established_years = sum(established_years),
    .groups = "drop"
  )

plot_dat <- combo_summary |>
  select(-any_of("n_runs")) |>
  left_join(denoms, by = c("R0", "K", "alpha", "seasonal_amp"))

## Floor for zeros: one event in the largest denominator, or half the
## smallest observed distance from 0 (or 1, for probabilities) if smaller
clamp_floor <- function(x, denom, prob = TRUE) {
  x <- x[is.finite(x)]
  gaps <- c(x[x > 0], if (prob) 1 - x[x < 1])
  min(1 / max(denom, na.rm = TRUE), min(gaps) / 2)
}

lo <- c(
  occ_median = clamp_floor(plot_dat$occ_median, settings$n_patch),
  p_persist = clamp_floor(plot_dat$p_persist, plot_dat$n_runs),
  p_burnout = clamp_floor(plot_dat$p_burnout, plot_dat$n_post_major),
  fadeout_hazard = clamp_floor(plot_dat$fadeout_hazard,
                               plot_dat$established_years, prob = FALSE)
)

plot_dat <- plot_dat |>
  mutate(
    across(c(occ_median, p_persist, p_burnout),
           \(x) pmin(pmax(x, lo[cur_column()]), 1 - lo[cur_column()])),
    fadeout_hazard = pmax(fadeout_hazard, lo["fadeout_hazard"])
  )

## White-to-red fill: logit scale for probabilities, log10 for rates
fill_scale <- function(var, fill_name) {
  prob <- var != "fadeout_hazard"
  if (prob) {
    limits <- c(lo[[var]], 1 - lo[[var]])
    breaks <- c(1e-4, 1e-3, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999, 1 - 1e-4)
  } else {
    limits <- c(lo[[var]], max(plot_dat[[var]], na.rm = TRUE))
    breaks <- 10^seq(-6, 2)
  }
  ## Drop interior breaks within 8% of the scale of a limit, so their labels
  ## do not collide with the limit labels
  tr <- scales::as.transform(if (prob) "logit" else "log10")
  tlim <- tr$transform(limits)
  gap <- 0.08 * diff(tlim)
  tb <- tr$transform(breaks)
  breaks <- breaks[tb > tlim[1] + gap & tb < tlim[2] - if (prob) gap else 0]
  ## The clamped limits stand for observed zeros (and ones)
  limit_labels <- function(b) {
    ifelse(abs(b - limits[1]) < 1e-12, "0",
           ifelse(prob & abs(b - limits[2]) < 1e-12, "1",
                  ## enough digits near 1 to show the distance from 1
                  vapply(b, \(v) format(
                    v,
                    digits = max(2, if (v < 1) ceiling(-log10(1 - v))),
                    drop0trailing = TRUE
                  ), character(1))))
  }
  scale_fill_distiller(
    name = fill_name, palette = "Reds", direction = 1,
    transform = tr,
    limits = limits, breaks = c(limits[1], breaks, if (prob) limits[2]),
    labels = limit_labels,
    na.value = "grey60"
  )
}

clamp_caption <- function(var) {
  if (var == "fadeout_hazard") {
    sprintf("zero drawn as %.2g; grey: no established episodes", lo[[var]])
  } else {
    sprintf("0 and 1 drawn as %.2g and 1 - %.2g; grey: undefined",
            lo[[var]], lo[[var]])
  }
}

heatmap <- function(var, fill_name, title) {
  ggplot(plot_dat,
         aes(x = factor(K / 1000), y = factor(R0), fill = .data[[var]])) +
    geom_raster() +
    facet_grid(seasonal_amp ~ alpha, labeller = facet_labels) +
    fill_scale(var, fill_name) +
    coord_cartesian(expand = FALSE) +
    labs(
      x = "Carrying capacity per patch (K, thousands)",
      y = expression(R[0]),
      title = title,
      subtitle = sprintf(
        "%d replicate(s) per cell, %d patches, %d years",
        settings$n_reps, settings$n_patch, settings$n_years
      ),
      caption = clamp_caption(var)
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
    data = filter(plot_dat, p_persist < 1 - lo["p_persist"]),
    shape = 4,
    colour = "black"
  ) +
  labs(caption = paste0(
    clamp_caption("occ_median"),
    "\nx: at least one replicate went globally extinct"
  ))

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
  "Fade-out hazard after establishment"
)

fig_width <- 16
fig_height <- 12

save_figure(p_occ, "heatmap_occupancy", fig_width, fig_height)
save_figure(p_persist, "heatmap_persistence", fig_width, fig_height)
save_figure(p_burnout, "heatmap_burnout", fig_width, fig_height)
save_figure(p_hazard, "heatmap_fadeout_hazard", fig_width, fig_height)


## ------------------------------------------------------------
## Amplitude (x) by R0 (y), faceted by K, at a single alpha
## ------------------------------------------------------------

## alpha value closest to 1e-4/day (on a log scale) in the grid
alpha_vals <- sort(unique(plot_dat$alpha))
amp_alpha <- alpha_vals[which.min(abs(log(alpha_vals / 1e-4)))]

amp_heatmap <- function(var, fill_name, title) {
  ggplot(filter(plot_dat, alpha == amp_alpha),
         aes(x = factor(seasonal_amp), y = factor(R0), fill = .data[[var]])) +
    geom_raster() +
    facet_wrap(~ K, nrow = 1,
               labeller = labeller(K = function(x) paste0("K = ", x))) +
    fill_scale(var, fill_name) +
    coord_cartesian(expand = FALSE) +
    labs(
      x = "Seasonal amplitude",
      y = expression(R[0]),
      title = title,
      subtitle = sprintf(
        "alpha = %g/day, %d replicate(s) per cell, %d patches, %d years",
        amp_alpha, settings$n_reps, settings$n_patch, settings$n_years
      ),
      caption = clamp_caption(var)
    ) +
    zmargin
}

amp_width <- 16
amp_height <- 4

save_figure(
  amp_heatmap("p_burnout", "P(burnout)", "Burnout probability vs. seasonality"),
  "heatmap_amp_burnout", amp_width, amp_height
)
save_figure(
  amp_heatmap("fadeout_hazard", "fade-outs per\npatch-year",
              "Fade-out hazard vs. seasonality"),
  "heatmap_amp_fadeout_hazard", amp_width, amp_height
)

cat("Output: ", outdir, "\n", sep = "")
