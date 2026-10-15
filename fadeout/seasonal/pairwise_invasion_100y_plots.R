## Figures for the pairwise invasion runs.
##
## Reads data/combo_summary.csv, data/run_summary.csv,
## data/occupancy_thinned.rds and data/settings.csv written by
## pairwise_invasion_100y.R and draws:
##
##   pip_outcomes      resident R0 (x) by invader R0 (y) heatmaps of outcome
##                     probabilities and late invader occupancy
##   invader_survival  fraction of runs with the invader still present, by
##                     years since invasion
##   occupancy_series  patch occupancy of each strain around the invasion
##
## Usage:
##   Rscript fadeout/seasonal/pairwise_invasion_100y_plots.R [--pilot]
##     [--susc-eq]
## --pilot and --susc-eq plot the output of runs made with those flags.

library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

theme_set(theme_bw())

zmargin <- theme(panel.spacing = grid::unit(0, "pt"))

pilot <- "--pilot" %in% commandArgs(trailingOnly = TRUE)

outdir <- here::here(
  "fadeout", "output",
  paste0(
    if (pilot) "pairwise_invasion_100y_pilot" else "pairwise_invasion_100y",
    if ("--susc-eq" %in% commandArgs(trailingOnly = TRUE)) "_susc_eq"
  )
)

dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)

combo_summary <- read.csv(file.path(outdir, "data", "combo_summary.csv"))
run_summary <- read.csv(file.path(outdir, "data", "run_summary.csv"))
occupancy_thinned <- readRDS(file.path(outdir, "data", "occupancy_thinned.rds"))

settings <- read.csv(file.path(outdir, "data", "settings.csv"))
settings <- setNames(as.list(settings$value), settings$name)

regime <- sprintf(
  "K = %g, alpha = %g/day, seasonal amplitude = %g; %g replicates per pair",
  settings$K, settings$alpha, settings$seasonal_amp, settings$n_reps
)

## Save each figure as both PDF and 300-dpi PNG (for HTML output)
save_figure <- function(plot, name, width, height) {
  for (ext in c("pdf", "png")) {
    ggsave(file.path(outdir, "figures", paste0(name, ".", ext)),
           plot, width = width, height = height, dpi = 300)
  }
}

## Resident and invader colours
strain_cols <- c(resident = "#2a78d6", invader = "#eb6834")

R0_res_lab <- function(x) paste0("resident R0 = ", x)
R0_inv_lab <- function(x) paste0("invader R0 = ", x)


## ------------------------------------------------------------
## Pairwise invasion plots
## ------------------------------------------------------------

pip_vars <- c(
  p_invader_persist = sprintf("P(invader present after %g y)",
                              settings$invasion_years),
  p_replace = "P(invader replaces resident)",
  p_coexist = "P(both strains present)",
  occ2_late = sprintf("invader occupancy, last %g y", settings$summary_years)
)

pip_dat <- combo_summary |>
  select(R0_res, R0_inv, all_of(names(pip_vars))) |>
  pivot_longer(all_of(names(pip_vars))) |>
  mutate(
    name = factor(pip_vars[name], levels = pip_vars),
    across(c(R0_res, R0_inv), factor)
  )

## Cells with no valid runs (resident always extinct before the invasion)
## are NaN and drawn grey without a label
any_absent <- any(is.na(pip_dat$value))

gg_pip <- ggplot(pip_dat, aes(R0_res, R0_inv, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(value), "", sprintf("%.2f", value)),
                colour = !is.na(value) & value > 0.6),
            size = 2.6, show.legend = FALSE) +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white")) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "grey40") +
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0, 1),
                       name = NULL, na.value = "grey80") +
  facet_wrap(~ name, nrow = 2) +
  coord_fixed(expand = FALSE) +
  labs(x = expression("Resident " * R[0]), y = expression("Invader " * R[0]),
       subtitle = regime,
       caption = if (any_absent)
         "Grey: resident globally extinct before the invasion in every run") +
  theme(panel.grid = element_blank())

save_figure(gg_pip, "pip_outcomes", width = 10, height = 9.5)


## ------------------------------------------------------------
## Invader survival after invasion
## ------------------------------------------------------------

## Kaplan-Meier is not needed: every run is followed for the full phase, so
## the curve is the empirical fraction still present. Runs where the
## resident was already extinct are dropped.
surv_times <- c(0, 10^seq(0, log10(settings$invasion_years * 365),
                          length.out = 200))

surv_dat <- run_summary |>
  filter(outcome != "resident_absent") |>
  select(R0_res, R0_inv, invader_extinction_day) |>
  cross_join(tibble(day = surv_times)) |>
  group_by(R0_res, R0_inv, day) |>
  summarise(
    surviving = mean(is.na(invader_extinction_day) |
                       invader_extinction_day > day),
    .groups = "drop"
  ) |>
  mutate(R0_inv = factor(R0_inv, levels = sort(unique(run_summary$R0_inv))))

## Ordered invader R0: one blue ramp, light to dark
inv_ramp <- colorRampPalette(c("#9ec5f4", "#3987e5", "#0d366b"))(
  nlevels(surv_dat$R0_inv)
)

gg_surv <- ggplot(surv_dat, aes(pmax(day, 1) / 365, surviving,
                                colour = R0_inv)) +
  geom_step(linewidth = 0.7) +
  scale_x_log10(labels = scales::label_number(drop0trailing = TRUE)) +
  scale_colour_manual(values = inv_ramp, name = expression("Invader " * R[0]),
                      drop = FALSE) +
  facet_wrap(~ R0_res, labeller = labeller(R0_res = R0_res_lab), nrow = 2) +
  labs(x = "Years since invasion",
       y = "Fraction of runs with the invader present",
       subtitle = regime) +
  zmargin

## Skipped when the resident was extinct before every invasion (some pilots)
if (nrow(surv_dat)) {
  save_figure(gg_surv, "invader_survival", width = 9, height = 6)
}


## ------------------------------------------------------------
## Occupancy around the invasion
## ------------------------------------------------------------

## Up to 20 replicates per pair
max_lines <- 20

series_dat <- occupancy_thinned |>
  inner_join(
    run_summary |>
      filter(replicate <= max_lines) |>
      select(run_id, R0_res, R0_inv, replicate),
    by = "run_id"
  ) |>
  transmute(
    run_id, R0_res, R0_inv,
    year = day / 365,
    resident = occ1 / settings$n_patch,
    invader = ifelse(day >= settings$resident_years * 365,
                     occ2 / settings$n_patch, NA)
  ) |>
  pivot_longer(c(resident, invader), names_to = "strain",
               values_to = "occupancy") |>
  filter(!is.na(occupancy)) |>
  mutate(strain = factor(strain, levels = names(strain_cols)))

gg_series <- ggplot(series_dat,
                    aes(year, occupancy, colour = strain,
                        group = interaction(run_id, strain))) +
  geom_vline(xintercept = settings$resident_years, colour = "grey50",
             linetype = 2) +
  ## Thousands of overlapping lines: rasterize them to keep the PDF small
  ggrastr::rasterise(geom_line(linewidth = 0.25, alpha = 0.5), dpi = 300) +
  scale_colour_manual(values = strain_cols, name = NULL) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1,
                                                   alpha = 1))) +
  facet_grid(R0_inv ~ R0_res,
             labeller = labeller(R0_res = R0_res_lab, R0_inv = R0_inv_lab)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = "Year", y = "Fraction of patches occupied",
       subtitle = sprintf("%s; %d replicates shown", regime,
                          min(max_lines, settings$n_reps))) +
  theme(legend.position = "top",
        strip.text.y = element_text(angle = 0)) +
  zmargin

## About 2.3 x 1.6 inches per panel
n_R0 <- n_distinct(series_dat$R0_res)
save_figure(gg_series, "occupancy_series",
            width = 2 + 2.3 * n_R0, height = 1.5 + 1.6 * n_R0)

cat("Figures: ", file.path(outdir, "figures"), "\n", sep = "")
