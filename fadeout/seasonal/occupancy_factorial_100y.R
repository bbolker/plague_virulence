## Seasonal vs constant-transmission patch occupancy over 100 years,
## factorial over R0 x K.
##
## Extends occupancy_seasonal_vs_constant_100y.R to a 3 x 3 grid of mean R0
## and carrying capacity K, run for both seasonal_amp = 0.40 and 0. Demography
## and initial conditions follow the 'ratchet' runs in single_strain_metapop.R
## (r = 0.125 per infectious period; each patch starts at a random point
## between a fresh outbreak and the endemic equilibrium), but time is in days
## with gamma = 0.2, dt = 1 day.
##
## Also crossed with three values of the between-patch transmission rate
## alpha. Each (R0, K, alpha, seasonal_amp) combination is a single
## realization; jobs are spread over a fixed number of forked workers.
##
## The state is thinned to every THIN-th day before conv_odin(), so one-step
## local-extinction and recolonization counts are meaningless and are not
## computed here.

library(plagueMetapop)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(odin)
library(future)
library(furrr)

theme_set(theme_bw())

zmargin <- theme(panel.spacing = grid::unit(0, "pt"))

## Okabe-Ito minus black and yellow
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

## 101 matches the seed of the ratchet runs; every job resets to it so that
## all combinations share the same uniform draws for their initial conditions
seed <- 101
set.seed(seed)

RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)


## ------------------------------------------------------------
## Parameters
## ------------------------------------------------------------

dt <- 1                   ## days
gamma <- 0.2              ## per day
r <- 0.125 * gamma        ## 0.125 per infectious period
n_patch <- 200
season_period <- 365
peak_day <- 15

## Infected count of a freshly invaded patch in the ratchet initial conditions
I_outbreak <- 10

n_years <- 100
t_max <- n_years * 365
nt <- round(t_max / dt)

## Keep every THIN-th day of output
THIN <- 10

## Amplitudes to compare; names become the scenario factor levels
amps <- c(seasonal = 0.40, constant = 0.00)

## Between-patch transmission, per day: 2e-5 is the ratchet runs' 1e-4 per
## infectious period, 1e-4 is the original 100-year run
alpha_vals <- c(2e-5, 5e-5, 1e-4)

## Jobs are spread over at most this many forked workers
n_workers <- 27

param_grid <- expand_grid(
  R0 = c(2, 2.5, 3),
  K = c(3e3, 1e4, 3e4),
  alpha = alpha_vals,
  scenario = names(amps)
) |>
  mutate(seasonal_amp = amps[scenario])


## ------------------------------------------------------------
## Compile the shared odin model
## ------------------------------------------------------------

model_file <- here::here("fadeout", "seasonal", "seasonal_model_metapop.R")

if (!file.exists(model_file)) {
  stop("Odin model file does not exist: ", model_file)
}

## odin::odin() substitutes its argument, so the path must be a variable.
## The compiled model is loaded into this process and reaches the workers by
## forking, which is why plan(multicore) is used rather than multisession.
gen <- suppressMessages(odin::odin(model_file))


## ------------------------------------------------------------
## One realization
## ------------------------------------------------------------

run_scenario <- function(R0, K, alpha, seasonal_amp) {

  ## Mersenne-Twister, as in the ratchet runs; furrr hands each worker an
  ## L'Ecuyer-CMRG stream, which this replaces
  set.seed(seed, kind = "Mersenne-Twister")

  beta0 <- R0 * gamma

  ## Deterministic equilibrium at the mean transmission rate
  S_star <- gamma * K / beta0
  I_star <- r * (K - S_star) / beta0

  ## Ratchet initial conditions: u = 1 is a fresh outbreak in a fully
  ## susceptible patch, u = 0 is the endemic equilibrium
  u <- runif(n_patch, min = 0, max = 1)
  S_ini <- round(K * u + S_star * (1 - u))
  I_ini <- cbind(round(I_outbreak * u + I_star * (1 - u)), rep(0, n_patch))

  mod <- gen$new(
    beta = c(beta0, 0),
    gamma = c(gamma, gamma),
    dt = dt,
    I_ini = I_ini,
    S_ini = S_ini,
    I2_ini = rep(0, n_patch),
    alpha = alpha,
    strain2_delay = .Machine$integer.max,
    r = rep(r, n_patch),
    K = rep(K, n_patch),
    season_period = season_period,
    seasonal_amp = seasonal_amp,
    peak_day = peak_day,
    n_patch = n_patch
  )

  t0 <- Sys.time()
  raw <- mod$run(seq(0L, nt))
  raw[, "step"] <- raw[, "step"] * dt

  ## Thin before conv_odin(), which is where the time goes
  keep <- seq(1, nrow(raw), by = THIN)
  runs <- conv_odin(raw[keep, , drop = FALSE])
  rm(raw)

  runs |>
    filter(state %in% c("S", "I1")) |>
    select(step, patch, state, value) |>
    pivot_wider(names_from = state, values_from = value) |>
    rename(I = I1) |>
    mutate(N = S + I, occupied = as.integer(I > 0)) |>
    group_by(step) |>
    summarise(
      occupied_patches = sum(occupied),
      global_I = sum(I),
      mean_N = mean(N),
      .groups = "drop"
    ) |>
    mutate(
      years = step / 365,
      elapsed_sec = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )
}


## ------------------------------------------------------------
## Run all combinations in parallel
## ------------------------------------------------------------

plan(multicore(workers = min(n_workers, nrow(param_grid))))

t_start <- Sys.time()

meta_summary <- param_grid |>
  mutate(
    res = future_pmap(
      list(R0, K, alpha, seasonal_amp),
      run_scenario,
      .options = furrr_options(seed = TRUE)
    )
  ) |>
  unnest(res) |>
  mutate(scenario = factor(scenario, levels = names(amps)))

plan(sequential)

cat(
  "total elapsed: ",
  round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1),
  " s\n",
  sep = ""
)


## ------------------------------------------------------------
## Summaries
## ------------------------------------------------------------

burnin <- 365

timing <- meta_summary |>
  distinct(R0, K, alpha, scenario, elapsed_sec)

print(timing, n = Inf)

overall_summary <- meta_summary |>
  filter(step >= burnin) |>
  group_by(R0, K, alpha, scenario) |>
  summarise(
    seasonal_amp = first(seasonal_amp),
    globally_persistent = all(global_I > 0),
    extinction_year =
      if (any(global_I == 0)) min(years[global_I == 0]) else NA_real_,
    min_occupied_patches = min(occupied_patches),
    mean_occupied_patches = mean(occupied_patches),
    max_occupied_patches = max(occupied_patches),
    mean_patch_population = mean(mean_N),
    .groups = "drop"
  )

print(overall_summary, n = Inf)

## Mean occupancy per decade, as a coarse trend check
decade_summary <- meta_summary |>
  filter(step >= burnin) |>
  mutate(decade = floor(years / 10) * 10) |>
  group_by(R0, K, alpha, scenario, decade) |>
  summarise(mean_occupied_patches = mean(occupied_patches), .groups = "drop")

print(
  decade_summary |>
    pivot_wider(names_from = decade, values_from = mean_occupied_patches),
  n = Inf,
  width = Inf
)


## ------------------------------------------------------------
## Outputs
## ------------------------------------------------------------

outdir <- here::here("fadeout", "output", "occupancy_factorial_100y")

dir.create(file.path(outdir, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)

write.csv(meta_summary, file.path(outdir, "data", "occupancy_thinned.csv"),
          row.names = FALSE)
write.csv(overall_summary, file.path(outdir, "data", "overall_summary.csv"),
          row.names = FALSE)
write.csv(decade_summary, file.path(outdir, "data", "decade_summary.csv"),
          row.names = FALSE)

## Zero spacing within each alpha group of K columns, a gap between groups
n_K <- n_distinct(param_grid$K)
alpha_gaps <- theme(
  panel.spacing.x = grid::unit(
    rep(c(rep(0, n_K - 1), 8), length.out = n_K * length(alpha_vals) - 1),
    "pt"
  )
)

## Save each figure as both PDF and 300-dpi PNG (for HTML output)
save_figure <- function(plot, name, width, height) {
  for (ext in c("pdf", "png")) {
    ggsave(file.path(outdir, "figures", paste0(name, ".", ext)),
           plot, width = width, height = height, dpi = 300)
  }
}

scenario_labels <- c(
  seasonal = sprintf("seasonal\n(amplitude = %g)", amps[["seasonal"]]),
  constant = "constant"
)

occupancy_ylab <- sprintf("Number of infected patches (n=%d)", n_patch)

facet_labels <- labeller(
  R0 = function(x) paste("R0 =", x),
  K = function(x) paste("K =", format(as.numeric(x), big.mark = ",")),
  alpha = function(x) paste0("alpha = ", x, "/day")
)


## ------------------------------------------------------------
## Plot 1: 100-year occupancy, R0 (rows) x alpha, K (columns)
## ------------------------------------------------------------

p_100y <- ggplot(
  meta_summary,
  aes(x = years, y = occupied_patches, colour = scenario)
) +
  geom_line(linewidth = 0.3) +
  facet_grid(R0 ~ alpha + K, labeller = facet_labels) +
  scale_colour_manual(values = okabe_ito, labels = scenario_labels, name = NULL) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1))) +
  scale_x_continuous(breaks = seq(0, n_years, by = 50)) +
  labs(
    x = "Time (years)",
    y = occupancy_ylab,
    title = sprintf(
      "Patch occupancy over %d years: seasonal vs constant transmission",
      n_years
    ),
    subtitle = sprintf(
      paste0(
        "r = %g/day, %d patches, ratchet initial conditions, ",
        "single realization per panel (seed %d), every %dth day plotted"
      ),
      r, n_patch, seed, THIN
    )
  ) +
  zmargin +
  alpha_gaps +
  theme(legend.position = "top")

save_figure(p_100y, "occupancy_factorial_100y", width = 20, height = 8)


## ------------------------------------------------------------
## Plot 2: first ten years, with the seasonal R0 curve
## ------------------------------------------------------------

zoom_years <- 10

zoom_occ <- meta_summary |>
  filter(years <= zoom_years)

seasonal_curve <- expand_grid(
  step = seq(0, zoom_years * 365, by = THIN),
  R0 = unique(param_grid$R0)
) |>
  mutate(
    R0_t = R0 * (1 + amps[["seasonal"]] *
                   cos(2 * pi * (step - peak_day) / season_period)),
    years = step / 365
  )

R0_colour <- "grey50"

## Secondary axis for R0, shared across panels
R0_scale <- n_patch / max(seasonal_curve$R0_t)

p_zoom <- ggplot(zoom_occ, aes(x = years)) +
  geom_line(
    data = seasonal_curve,
    aes(y = R0_t * R0_scale),
    colour = R0_colour,
    linewidth = 0.3,
    show.legend = FALSE
  ) +
  geom_line(aes(y = occupied_patches, colour = scenario), linewidth = 0.3) +
  facet_grid(R0 ~ alpha + K, labeller = facet_labels) +
  scale_colour_manual(values = okabe_ito, labels = scenario_labels, name = NULL) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1))) +
  scale_y_continuous(
    name = occupancy_ylab,
    sec.axis = sec_axis(~ . / R0_scale, name = "Seasonal R0 (grey)")
  ) +
  scale_x_continuous(breaks = seq(0, zoom_years, by = 5)) +
  labs(
    x = "Time (years)",
    title = sprintf(
      "First %d years: occupancy and seasonal R0",
      zoom_years
    )
  ) +
  zmargin +
  alpha_gaps +
  theme(
    legend.position = "top",
    axis.title.y.right = element_text(colour = R0_colour),
    axis.text.y.right = element_text(colour = R0_colour)
  )

save_figure(p_zoom, "occupancy_factorial_first10y", width = 20, height = 8)


## ------------------------------------------------------------
## Plot 3: occupancy distribution over the last 50 years of each run
## ------------------------------------------------------------

summary_years <- 50

late_summary <- meta_summary |>
  filter(years >= n_years - summary_years) |>
  group_by(R0, K, alpha, scenario) |>
  summarise(
    median = median(occupied_patches),
    lwr = quantile(occupied_patches, 0.025),
    upr = quantile(occupied_patches, 0.975),
    .groups = "drop"
  )

write.csv(late_summary, file.path(outdir, "data", "late_summary.csv"),
          row.names = FALSE)

p_late <- ggplot(
  late_summary,
  aes(
    x = factor(R0),
    y = median,
    ymin = lwr,
    ymax = upr,
    colour = scenario,
    shape = factor(alpha),
    group = interaction(scenario, alpha)
  )
) +
  geom_pointrange(position = position_dodge(width = 0.7), size = 0.3) +
  facet_grid(. ~ K, labeller = facet_labels) +
  scale_colour_manual(values = okabe_ito, labels = scenario_labels, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = "alpha (/day)") +
  labs(
    x = "R0",
    y = occupancy_ylab,
    title = sprintf(
      "Occupancy over the last %d years: median and 95%% range",
      summary_years
    )
  ) +
  zmargin +
  theme(legend.position = "top")

save_figure(p_late, "occupancy_factorial_late_summary", width = 8, height = 4)

cat("Output: ", outdir, "\n", sep = "")
