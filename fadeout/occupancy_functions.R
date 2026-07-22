## Reusable helpers for exploratory single-strain occupancy simulations.
##
## These functions deliberately preserve the explicit initialization used in
## fadeout/single_strain_metapop.R rather than using discrete_run() defaults.

make_fadeout_initial_conditions <- function(R0, K, r, gamma, n_patch,
                                            I_outbreak = 10, seed = 101,
                                            initialization = "interpolated") {
  set.seed(seed)

  eq <- plagueMetapop::ode_eq(
    beta = R0,
    gamma = gamma,
    K = K,
    r = r,
    logistic_growth = 1
  )

  S_star <- unname(eq["eq_S"])
  I_star <- unname(eq["eq_I"])
  initialization <- match.arg(
    initialization,
    c("interpolated", "virgin_soil")
  )

  if (initialization == "interpolated") {
    u <- stats::runif(n_patch, min = 0, max = 1)
    S_ini <- round(K * u + S_star * (1 - u))
    I1_ini <- round(I_outbreak * u + I_star * (1 - u))
  } else {
    I1_ini <- rep(as.integer(I_outbreak), n_patch)
    S_ini <- rep(as.integer(K - I_outbreak), n_patch)
  }

  list(
    S_ini = S_ini,
    I_ini = cbind(I1_ini, rep(0, n_patch)),
    S_star = S_star,
    I_star = I_star
  )
}

summarize_fadeout_occupancy <- function(runs, n_patch) {
  infected <- runs |>
    dplyr::filter(state == "I1") |>
    dplyr::select(step, patch, I = value) |>
    dplyr::arrange(patch, step) |>
    dplyr::mutate(occupied = as.integer(I > 0)) |>
    dplyr::group_by(patch) |>
    dplyr::arrange(step, .by_group = TRUE) |>
    dplyr::mutate(
      prev_occupied = dplyr::lag(occupied),
      local_extinction = as.integer(prev_occupied == 1 & occupied == 0),
      recolonization = as.integer(prev_occupied == 0 & occupied == 1)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      local_extinction = tidyr::replace_na(local_extinction, 0L),
      recolonization = tidyr::replace_na(recolonization, 0L)
    )

  meta_summary <- infected |>
    dplyr::group_by(step) |>
    dplyr::summarise(
      occupied_patches = sum(occupied),
      occupancy_fraction = occupied_patches / n_patch,
      global_I = sum(I),
      local_extinctions = sum(local_extinction),
      recolonizations = sum(recolonization),
      .groups = "drop"
    )

  list(infected = infected, meta_summary = meta_summary)
}

add_established_occupancy <- function(infected, tau, dt) {
  required <- c("step", "patch", "I")
  missing <- setdiff(required, names(infected))
  if (length(missing) > 0L) {
    stop("Patch trajectories are missing columns: ",
         paste(missing, collapse = ", "))
  }
  if (length(tau) != 1L || !is.finite(tau) || tau < 0) {
    stop("tau must be one finite, non-negative value")
  }
  if (length(dt) != 1L || !is.finite(dt) || dt <= 0) {
    stop("dt must be one finite, positive value")
  }

  window_steps <- as.integer(round(tau / dt))
  if (!isTRUE(all.equal(window_steps * dt, tau, tolerance = 1e-8))) {
    stop("tau must be an integer multiple of dt; got tau=", tau,
         " and dt=", dt)
  }

  infected |>
    dplyr::select(dplyr::all_of(required)) |>
    dplyr::arrange(patch, step) |>
    dplyr::group_by(patch) |>
    dplyr::group_modify(function(.x, .y) {
      n <- nrow(.x)
      established <- rep(NA_integer_, n)
      valid <- seq_len(max(0L, n - window_steps))
      if (length(valid) > 0L) {
        zero_cumulative <- c(0L, cumsum(.x$I <= 0))
        window_end <- valid + window_steps
        zero_count <- zero_cumulative[window_end + 1L] -
          zero_cumulative[valid]
        established[valid] <- as.integer(zero_count == 0L)
      }
      .x$occupied_raw <- as.integer(.x$I > 0)
      .x$occupied_established <- established
      .x
    }) |>
    dplyr::ungroup() |>
    dplyr::mutate(tau = tau, window_steps = window_steps)
}

summarize_established_occupancy <- function(infected, n_patch, tau, dt) {
  if (length(n_patch) != 1L || n_patch <= 0) {
    stop("n_patch must be one positive value")
  }

  established <- add_established_occupancy(infected, tau = tau, dt = dt)
  established |>
    dplyr::group_by(step, tau) |>
    dplyr::summarise(
      occupied_patches_raw = sum(occupied_raw),
      occupancy_fraction_raw = occupied_patches_raw / n_patch,
      occupied_patches_established = if (all(is.na(occupied_established))) {
        NA_integer_
      } else {
        sum(occupied_established, na.rm = TRUE)
      },
      patches_with_complete_window = sum(!is.na(occupied_established)),
      occupancy_fraction_established = if (
        patches_with_complete_window == 0L
      ) {
        NA_real_
      } else {
        occupied_patches_established / patches_with_complete_window
      },
      global_I = sum(I),
      .groups = "drop"
    )
}

classify_infection_episodes_with_censoring <- function(infected, tau, dt) {
  ## Corrected audit classification. Episodes are independent within a patch.
  ## A qualifying episode is persistent from its beginning; a completed
  ## non-qualifying episode is transient; an episode still infected at the
  ## simulation boundary without a complete qualifying window is censored.
  classified <- add_established_occupancy(infected, tau = tau, dt = dt) |>
    dplyr::group_by(patch) |>
    dplyr::arrange(step, .by_group = TRUE) |>
    dplyr::mutate(
      episode_start = occupied_raw == 1L &
        dplyr::lag(occupied_raw, default = 0L) == 0L,
      episode_id = cumsum(episode_start)
    ) |>
    dplyr::ungroup()

  simulation_end <- max(classified$step)
  episode_types <- classified |>
    dplyr::filter(occupied_raw == 1L) |>
    dplyr::group_by(patch, episode_id) |>
    dplyr::summarise(
      episode_persistent = any(occupied_established == 1L, na.rm = TRUE),
      episode_start_time = min(step),
      episode_end_time = max(step),
      episode_right_censored = episode_end_time == simulation_end,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      episode_type = dplyr::case_when(
        episode_persistent ~ "persistent",
        episode_right_censored ~ "censored",
        TRUE ~ "transient"
      )
    )

  classified |>
    dplyr::left_join(episode_types, by = c("patch", "episode_id")) |>
    dplyr::mutate(
      episode_type = dplyr::if_else(
        occupied_raw == 0L, "uninfected", episode_type
      )
    )
}

summarize_transient_source_pressure <- function(infected, n_patch, tau, dt) {
  classified <- classify_infection_episodes_with_censoring(
    infected, tau = tau, dt = dt
  )

  classified |>
    dplyr::group_by(step) |>
    dplyr::summarise(
      transient_I = sum(I[episode_type == "transient"]),
      persistent_I = sum(I[episode_type == "persistent"]),
      censored_I = sum(I[episode_type == "censored"]),
      classifiable_I = transient_I + persistent_I,
      total_I = classifiable_I + censored_I,
      transient_source_share = dplyr::if_else(
        classifiable_I > 0, transient_I / classifiable_I, NA_real_
      ),
      persistent_source_share = dplyr::if_else(
        classifiable_I > 0, persistent_I / classifiable_I, NA_real_
      ),
      transient_share_of_total_I = dplyr::if_else(
        total_I > 0, transient_I / total_I, NA_real_
      ),
      censored_share_of_total_I = dplyr::if_else(
        total_I > 0, censored_I / total_I, NA_real_
      ),
      infected_patches = sum(occupied_raw),
      transient_patches = sum(episode_type == "transient"),
      persistent_patches = sum(episode_type == "persistent"),
      censored_patches = sum(episode_type == "censored"),
      transient_patch_fraction = dplyr::if_else(
        infected_patches > 0, transient_patches / infected_patches, NA_real_
      ),
      persistent_patch_fraction = dplyr::if_else(
        infected_patches > 0, persistent_patches / infected_patches, NA_real_
      ),
      censored_patch_fraction = dplyr::if_else(
        infected_patches > 0, censored_patches / infected_patches, NA_real_
      ),
      persistent_occupancy_fraction = persistent_patches / n_patch,
      .groups = "drop"
    )
}

summarize_occupancy_curve <- function(meta_summary, early_window = 100) {
  global_min_index <- which.min(meta_summary$occupied_patches)
  early_rows <- which(meta_summary$step <= early_window)
  early_min_index <- early_rows[
    which.min(meta_summary$occupied_patches[early_rows])
  ]
  after_early_min <- meta_summary[
    early_min_index:nrow(meta_summary), , drop = FALSE
  ]

  data.frame(
    initially_occupied = meta_summary$occupied_patches[1],
    early_min_occupied_patches =
      meta_summary$occupied_patches[early_min_index],
    early_min_occupancy_fraction =
      meta_summary$occupancy_fraction[early_min_index],
    time_of_early_minimum = meta_summary$step[early_min_index],
    global_min_occupied_patches =
      meta_summary$occupied_patches[global_min_index],
    global_min_occupancy_fraction =
      meta_summary$occupancy_fraction[global_min_index],
    time_of_global_minimum = meta_summary$step[global_min_index],
    final_occupied_patches = meta_summary$occupied_patches[nrow(meta_summary)],
    recovery_from_early_minimum =
      meta_summary$occupied_patches[nrow(meta_summary)] -
      meta_summary$occupied_patches[early_min_index],
    post_early_minimum_recolonizations =
      sum(after_early_min$recolonizations),
    globally_persistent = all(meta_summary$global_I > 0)
  )
}

apply_full_occupancy_absorbing_boundary <- function(meta_summary) {
  ## Full occupancy is treated as absorbing only after the curve has first
  ## fallen below full occupancy. This prevents step zero from triggering the
  ## boundary before the initial epidemic trough has occurred.
  n_patch <- max(meta_summary$occupied_patches)
  below_full <- which(meta_summary$occupied_patches < n_patch)
  absorption_step <- NA_real_

  if (length(below_full) > 0) {
    returned <- which(
      seq_len(nrow(meta_summary)) > below_full[1] &
        meta_summary$occupied_patches == n_patch
    )
    if (length(returned) > 0) {
      absorption_index <- returned[1]
      absorption_step <- meta_summary$step[absorption_index]
      meta_summary$occupied_patches[absorption_index:nrow(meta_summary)] <-
        n_patch
      meta_summary$occupancy_fraction[absorption_index:nrow(meta_summary)] <-
        1
    }
  }

  attr(meta_summary, "absorption_step") <- absorption_step
  meta_summary
}

run_fadeout_occupancy <- function(params, gen = NULL, keep_patch_data = FALSE) {
  required <- c(
    "R0", "K", "r", "alpha", "gamma", "n_patch", "dt", "t_max",
    "I_outbreak", "initialization_seed", "simulation_seed"
  )
  missing <- setdiff(required, names(params))
  if (length(missing) > 0) {
    stop("Missing parameters: ", paste(missing, collapse = ", "))
  }

  initial <- make_fadeout_initial_conditions(
    R0 = params$R0,
    K = params$K,
    r = params$r,
    gamma = params$gamma,
    n_patch = params$n_patch,
    I_outbreak = params$I_outbreak,
    seed = params$initialization_seed,
    initialization = if (is.null(params$initialization)) {
      "interpolated"
    } else {
      params$initialization
    }
  )

  if (is.null(gen)) {
    gen <- plagueMetapop::compile_odin("euler_odin_def.R")
  }

  set.seed(params$simulation_seed)
  sim <- gen$new(
    beta = c(params$R0, 0),
    gamma = c(params$gamma, params$gamma),
    dt = params$dt,
    r = rep(params$r, params$n_patch),
    K = rep(params$K, params$n_patch),
    S_ini = initial$S_ini,
    I_ini = initial$I_ini,
    I2_ini = rep(0, params$n_patch),
    alpha = params$alpha,
    n_patch = params$n_patch,
    strain2_delay = as.integer(.Machine$integer.max),
    logistic_growth = 1,
    reedfrost = 0
  )

  raw_runs <- sim$run(seq.int(0L, round(params$t_max / params$dt)))
  raw_runs[, "step"] <- raw_runs[, "step"] * params$dt
  runs <- plagueMetapop::conv_odin(raw_runs)
  occupancy <- summarize_fadeout_occupancy(runs, params$n_patch)

  list(
    params = params,
    equilibrium = c(S_star = initial$S_star, I_star = initial$I_star),
    initial = initial[c("S_ini", "I_ini")],
    infected = if (keep_patch_data) occupancy$infected else NULL,
    meta_summary = occupancy$meta_summary,
    curve_summary = summarize_occupancy_curve(occupancy$meta_summary)
  )
}
