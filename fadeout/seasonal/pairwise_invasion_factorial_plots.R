## Figures for the factorial pairwise invasion runs on Nibi.
##
## Reads the combined output of odin/sharcnet/pairwise_invasion_factorial_combine.R
## (odin/sharcnet/pairwise_invasion_factorial[_mini].rds) and writes to
## fadeout/output/pairwise_invasion_factorial[_mini]/:
##
##   data/combo_summary.csv
##   figures/pip_<metric>_amp<amplitude>   resident R0 (x) by invader R0 (y)
##       heatmaps, faceted alpha ~ K, one figure per metric and amplitude
##   figures/occupancy_series_K<K>_alpha<alpha>_amp<amplitude>   with --series
##
## Usage:
##   Rscript fadeout/seasonal/pairwise_invasion_factorial_plots.R [--mini]
##     [--series K,alpha,amplitude]
## e.g. --series 5000,5e-5,0.2

library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(optparse)

theme_set(theme_bw())

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-m", "--mini"), action = "store_true", default = FALSE,
              help = "plot the mini run"),
  make_option("--series", type = "character", default = NULL,
              help = "also draw the occupancy series for K,alpha,amplitude")
)))

base_fn <- "pairwise_factorial"
if (opt$mini) base_fn <- paste0(base_fn, "_mini")

res <- readRDS(here::here("odin", "sharcnet", paste0(base_fn, ".rds")))
settings <- res$settings

outdir <- here::here(
  "fadeout", "output",
  if (opt$mini) "pairwise_invasion_factorial_mini" else "pairwise_invasion_factorial"
)
dir.create(file.path(outdir, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)

write.csv(res$combo_summary, file.path(outdir, "data", "combo_summary.csv"),
          row.names = FALSE)

if (length(res$missing_tasks)) {
  cat("note:", length(res$missing_tasks), "tasks missing from the combined output\n")
}

## Save each figure as both PDF and PNG (for HTML output)
save_figure <- function(plot, name, width, height, dpi = 300) {
  for (ext in c("pdf", "png")) {
    ggsave(file.path(outdir, "figures", paste0(name, ".", ext)),
           plot, width = width, height = height, dpi = dpi)
  }
}

## Resolution of the large occupancy-series figure (rasterized lines and PNG)
series_dpi <- 200

facet_labels <- labeller(
  K = \(x) paste0("K = ", x),
  alpha = \(x) paste0("alpha = ", x, "/day")
)

n_reps <- max(res$combo_summary$n_runs)


## ------------------------------------------------------------
## Pairwise invasion plots: one figure per metric and amplitude
## ------------------------------------------------------------

pip_vars <- c(
  p_invader_persist = sprintf("P(invader present after %g y)",
                              settings$invasion_years),
  p_replace = "P(invader replaces resident)",
  p_coexist = "P(both strains present)",
  occ2_late = sprintf("Invader occupancy, last %g y", settings$summary_years)
)

pip_dat <- res$combo_summary |>
  mutate(across(c(R0_res, R0_inv), factor))

for (var in names(pip_vars)) {
  for (amp in sort(unique(pip_dat$seasonal_amp))) {
    d <- pip_dat |>
      filter(seasonal_amp == amp) |>
      mutate(value = .data[[var]])
    gg <- ggplot(d, aes(R0_res, R0_inv, fill = value)) +
      geom_tile(colour = "white", linewidth = 0.4) +
      geom_text(aes(label = ifelse(is.na(value), "", sprintf("%.2f", value)),
                    colour = !is.na(value) & value > 0.6),
                size = 1.9, show.legend = FALSE) +
      scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white")) +
      geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "grey40") +
      scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0, 1),
                           name = NULL, na.value = "grey80") +
      facet_grid(alpha ~ K, labeller = facet_labels) +
      coord_fixed(expand = FALSE) +
      labs(x = expression("Resident " * R[0]),
           y = expression("Invader " * R[0]),
           title = pip_vars[[var]],
           subtitle = sprintf(
             "Seasonal amplitude %g; %d replicates per pair; grey: resident extinct before every invasion",
             amp, n_reps)) +
      theme(panel.grid = element_blank())
    save_figure(gg, sprintf("pip_%s_amp%g", var, amp), width = 11, height = 10)
  }
}


## ------------------------------------------------------------
## Occupancy series for one environment (--series)
## ------------------------------------------------------------

if (!is.null(opt$series)) {
  env <- as.numeric(strsplit(opt$series, ",")[[1]])
  stopifnot(length(env) == 3)
  near <- \(x, y) abs(x - y) < 1e-9 * pmax(1, abs(y))

  ## task ids of this environment, with their R0 pair
  tasks <- res$run_summary |>
    filter(near(K, env[1]), near(alpha, env[2]), near(seasonal_amp, env[3])) |>
    distinct(task_id, R0_res, R0_inv)
  if (!nrow(tasks)) stop("no runs for K, alpha, amplitude = ", opt$series)

  strain_cols <- c(resident = "#2a78d6", invader = "#eb6834")
  series_dat <- res$thinned |>
    inner_join(tasks, by = "task_id") |>
    transmute(
      run = interaction(task_id, replicate),
      R0_res, R0_inv,
      year = day / 365,
      resident = occ1 / settings$n_patch,
      invader = ifelse(day >= settings$resident_years * 365,
                       occ2 / settings$n_patch, NA)
    ) |>
    pivot_longer(c(resident, invader), names_to = "strain",
                 values_to = "occupancy") |>
    filter(!is.na(occupancy)) |>
    mutate(strain = factor(strain, levels = names(strain_cols)))

  gg <- ggplot(series_dat,
               aes(year, occupancy, colour = strain,
                   group = interaction(run, strain))) +
    geom_vline(xintercept = settings$resident_years, colour = "grey50",
               linetype = 2) +
    ## Thousands of overlapping lines: rasterize them to keep the PDF small
    ggrastr::rasterise(geom_line(linewidth = 0.25, alpha = 0.5),
                       dpi = series_dpi) +
    scale_colour_manual(values = strain_cols, name = NULL) +
    guides(colour = guide_legend(override.aes = list(linewidth = 1,
                                                     alpha = 1))) +
    facet_grid(R0_inv ~ R0_res,
               labeller = labeller(R0_res = \(x) paste0("resident R0 = ", x),
                                   R0_inv = \(x) paste0("invader R0 = ", x))) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    labs(x = "Year", y = "Fraction of patches occupied",
         subtitle = sprintf("K = %g, alpha = %g/day, seasonal amplitude %g; saved replicates",
                            env[1], env[2], env[3])) +
    theme(legend.position = "top", strip.text.y = element_text(angle = 0),
          panel.spacing = grid::unit(0, "pt"))

  n_R0 <- n_distinct(tasks$R0_res)
  save_figure(gg, sprintf("occupancy_series_K%g_alpha%g_amp%g",
                          env[1], env[2], env[3]),
              width = 2 + 2.3 * n_R0, height = 1.5 + 1.6 * n_R0,
              dpi = series_dpi)
}

cat("Figures: ", file.path(outdir, "figures"), "\n", sep = "")
