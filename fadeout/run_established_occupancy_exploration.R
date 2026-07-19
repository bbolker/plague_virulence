library(plagueMetapop)
library(dplyr)
library(here)

source(here::here("fadeout", "occupancy_exploration_functions.R"))

## This deliberately duplicates the scenario specification in
## run_occupancy_exploration.R so changes here cannot overwrite its outputs.
baseline <- list(
  R0 = 2.5, K = 1e4, r = 0.125, alpha = 1e-4, gamma = 1,
  n_patch = 200, dt = 0.1, t_max = 2000, I_outbreak = 10,
  initialization = "virgin_soil",
  initialization_seed = 101, simulation_seed = 102
)
parameter_values <- list(
  R0 = c(1.5, 2, 2.5, 3, 4, 5),
  K = c(1000, 3000, 10000, 30000, 100000),
  r = c(0.05, 0.1, 0.125, 0.2, 0.4),
  alpha = c(1e-5, 3e-5, 1e-4, 3e-4, 1e-3)
)
## A single conservative window is used to avoid treating long transient
## infections as established. With dt=0.1, tau=50 spans 500 stored steps.
tau_values <- 50

make_scenarios <- function(baseline, parameter_values) {
  unlist(lapply(names(parameter_values), function(parameter) {
    lapply(parameter_values[[parameter]], function(value) {
      params <- baseline
      params[[parameter]] <- value
      list(varied_parameter = parameter, parameter_value = value, params = params)
    })
  }), recursive = FALSE)
}
last_nonmissing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_real_ else tail(x, 1L)
}
time_of_minimum <- function(time, value) {
  ok <- !is.na(value)
  if (!any(ok)) NA_real_ else time[ok][which.min(value[ok])]
}

outdir <- here::here("fadeout", "output", "patch_occupancy_established")
datadir <- file.path(outdir, "data")
dir.create(datadir, recursive = TRUE, showWarnings = FALSE)

scenarios <- make_scenarios(baseline, parameter_values)
gen <- compile_odin("euler_odin_def.R")
results <- vector("list", length(scenarios))

for (i in seq_along(scenarios)) {
  scenario <- scenarios[[i]]
  message(sprintf("Running established occupancy: %s = %g (%d/%d)",
                  scenario$varied_parameter, scenario$parameter_value,
                  i, length(scenarios)))
  simulation <- run_fadeout_occupancy(
    scenario$params, gen = gen, keep_patch_data = TRUE
  )
  summaries <- bind_rows(lapply(tau_values, function(tau) {
    summarize_established_occupancy(
      simulation$infected, scenario$params$n_patch, tau, scenario$params$dt
    )
  }))
  results[[i]] <- list(
    varied_parameter = scenario$varied_parameter,
    parameter_value = scenario$parameter_value,
    params = scenario$params,
    occupancy_summary = summaries
  )
  rm(simulation, summaries)
  invisible(gc())
}

curves <- bind_rows(lapply(results, function(x) {
  x$occupancy_summary |>
    mutate(varied_parameter = x$varied_parameter,
           parameter_value = x$parameter_value,
           R0 = x$params$R0, K = x$params$K, r = x$params$r,
           alpha = x$params$alpha, .before = 1)
}))

diagnostics <- curves |>
  group_by(varied_parameter, parameter_value, R0, K, r, alpha, tau) |>
  summarise(
    minimum_raw_occupancy = min(occupancy_fraction_raw),
    minimum_established_occupancy = min(occupancy_fraction_established, na.rm = TRUE),
    final_raw_occupancy = tail(occupancy_fraction_raw, 1L),
    final_established_occupancy = last_nonmissing(occupancy_fraction_established),
    mean_difference_raw_minus_established = mean(
      occupancy_fraction_raw - occupancy_fraction_established, na.rm = TRUE),
    max_difference_raw_minus_established = max(
      occupancy_fraction_raw - occupancy_fraction_established, na.rm = TRUE),
    time_of_minimum_raw = step[which.min(occupancy_fraction_raw)],
    time_of_minimum_established = time_of_minimum(step, occupancy_fraction_established),
    number_of_nonmissing_time_points = sum(!is.na(occupancy_fraction_established)),
    .groups = "drop"
  )

saveRDS(results, file.path(datadir, "established_occupancy_results.rds"))
write.csv(curves, file.path(datadir, "established_occupancy_curves.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(datadir, "established_occupancy_diagnostics.csv"), row.names = FALSE)
message("Finished. New outputs saved in: ", outdir)
