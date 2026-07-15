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
