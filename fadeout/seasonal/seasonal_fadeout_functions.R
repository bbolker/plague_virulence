## Reusable diagnostics for recurrent fade-out in the seasonal,
## single-strain stochastic metapopulation model.

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Required package(s) not installed: ", paste(missing, collapse = ", "))
  }
}

calculate_intrinsic_period <- function(R0, gamma, r, tolerance = 1e-10) {
  values <- c(R0 = R0, gamma = gamma, r = r)
  if (any(!is.finite(values)) || gamma <= 0 || r <= 0) {
    stop("R0, gamma, and r must be finite; gamma and r must be positive")
  }
  if (R0 <= 1) stop("Intrinsic period requires R0 > 1")

  omega_squared <- gamma * r * (R0 - 1) / R0 -
    r^2 / (4 * R0^2)
  if (omega_squared <= 0) {
    stop("Intrinsic-period expression is not positive: ", omega_squared)
  }

  jacobian <- matrix(
    c(-r / R0, r * (R0 - 1) / R0, -gamma, 0),
    nrow = 2,
    byrow = FALSE
  )
  eigenvalues <- eigen(jacobian, only.values = TRUE)$values
  eigen_omega <- max(abs(Im(eigenvalues)))
  analytic_omega <- sqrt(omega_squared)
  if (!is.finite(eigen_omega) ||
      abs(eigen_omega - analytic_omega) >
        tolerance * max(1, analytic_omega)) {
    stop(
      "Analytic and Jacobian intrinsic frequencies disagree: ",
      analytic_omega, " versus ", eigen_omega
    )
  }

  list(
    T0 = 2 * pi / analytic_omega,
    omega0 = analytic_omega,
    eigen_omega0 = eigen_omega,
    eigenvalues = eigenvalues,
    jacobian = jacobian
  )
}

validate_seasonal_grid <- function(grid) {
  required <- c(
    "combo_id", "R0", "gamma", "r", "K", "alpha", "seasonal_amp",
    "season_period", "peak_day", "n_patch", "dt", "t_max", "n_reps",
    "base_seed", "threshold_multiplier"
  )
  missing <- setdiff(required, names(grid))
  if (length(missing)) {
    stop("Parameter grid is missing: ", paste(missing, collapse = ", "))
  }
  if (!nrow(grid)) stop("Parameter grid has no rows")
  if (anyDuplicated(grid$combo_id)) stop("combo_id values must be unique")

  numeric_columns <- setdiff(required, "combo_id")
  for (column in numeric_columns) {
    grid[[column]] <- suppressWarnings(as.numeric(grid[[column]]))
    if (any(!is.finite(grid[[column]]))) {
      stop("Column ", column, " contains non-numeric or non-finite values")
    }
  }
  if (any(grid$R0 <= 1)) stop("All R0 values must be > 1")
  if (any(grid$gamma <= 0 | grid$r <= 0 | grid$K <= 0 |
          grid$alpha < 0 | grid$season_period <= 0 | grid$dt <= 0 |
          grid$t_max <= 0 | grid$n_patch < 1 | grid$n_reps < 1 |
          grid$threshold_multiplier <= 0)) {
    stop("Grid contains an invalid non-positive parameter")
  }
  if (any(grid$seasonal_amp < 0 | grid$seasonal_amp > 1)) {
    stop("seasonal_amp must lie in [0, 1] so beta(t) cannot be negative")
  }
  integer_columns <- c("n_patch", "n_reps", "base_seed")
  for (column in integer_columns) {
    if (any(grid[[column]] != round(grid[[column]]))) {
      stop(column, " must contain integers")
    }
    grid[[column]] <- as.integer(grid[[column]])
  }
  grid
}

make_parameter_label <- function(row) {
  sprintf(
    "%s | R0=%g, gamma=%g, r=%g, K=%g, alpha=%g, amp=%g, period=%g",
    row$combo_id, row$R0, row$gamma, row$r, row$K, row$alpha,
    row$seasonal_amp, row$season_period
  )
}

safe_file_id <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
}

plot_occupancy_raster <- function(
    patch_state, parameters, replicate, seed, T0, fadeout_threshold,
    figure_dir) {
  raster_data <- patch_state |>
    dplyr::transmute(
      time = .data$time,
      patch = factor(.data$patch),
      occupied = factor(
        as.integer(.data$I > 0),
        levels = c(0, 1),
        labels = c("Uninfected", "Infected")
      )
    )
  plot <- ggplot2::ggplot(
    raster_data,
    ggplot2::aes(.data$time, .data$patch, fill = .data$occupied)
  ) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_manual(
      values = c("Uninfected" = "white", "Infected" = "black"),
      breaks = "Infected",
      labels = "Infected (I > 0)",
      drop = FALSE
    ) +
    ggplot2::labs(
      x = "Time (days)", y = "Patch", fill = NULL,
      title = "Seasonal single-strain patch occupancy",
      subtitle = paste0(
        make_parameter_label(parameters),
        sprintf(
          "\nreplicate=%d; seed=%d; T0=%.2f days; fade-out threshold=%.2f days",
          replicate, seed, T0, fadeout_threshold
        )
      )
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      legend.position = "top",
      plot.subtitle = ggplot2::element_text(size = 8)
    )
  if (parameters$n_patch > 40) {
    plot <- plot + ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  }
  filename <- sprintf(
    "occupancy_raster_%s_rep%03d.pdf",
    safe_file_id(parameters$combo_id), replicate
  )
  ggplot2::ggsave(
    file.path(figure_dir, filename), plot,
    width = 10, height = if (parameters$n_patch <= 20) 5 else 7
  )
  invisible(filename)
}

extract_patch_state <- function(runs) {
  required <- c("step", "patch", "state", "value")
  if (!all(required %in% names(runs))) {
    stop("Converted odin trajectory lacks: ",
         paste(setdiff(required, names(runs)), collapse = ", "))
  }
  state <- runs |>
    dplyr::filter(.data$state %in% c("S", "I1")) |>
    dplyr::select(
      time = "step", patch = "patch",
      state = "state", value = "value"
    ) |>
    tidyr::pivot_wider(names_from = "state", values_from = "value") |>
    dplyr::rename(I = "I1") |>
    dplyr::arrange(.data$patch, .data$time) |>
    dplyr::mutate(occupied = .data$I > 0)
  if (anyNA(state$I)) stop("Missing infected counts after reshaping trajectory")
  state
}

extract_infection_episodes <- function(
    patch_state, combo_id, replicate, T0, threshold_multiplier) {
  required <- c("time", "patch", "I")
  if (!all(required %in% names(patch_state))) {
    stop("patch_state lacks: ",
         paste(setdiff(required, names(patch_state)), collapse = ", "))
  }
  threshold <- threshold_multiplier * T0
  censoring_time <- max(patch_state$time)

  per_patch <- split(patch_state, patch_state$patch)
  episodes <- lapply(per_patch, function(dat) {
    dat <- dat[order(dat$time), , drop = FALSE]
    occupied <- dat$I > 0
    runs <- rle(occupied)
    run_end <- cumsum(runs$lengths)
    run_start <- run_end - runs$lengths + 1L
    infected_runs <- which(runs$values)
    if (!length(infected_runs)) return(NULL)

    do.call(rbind, lapply(seq_along(infected_runs), function(episode_id) {
      run_index <- infected_runs[episode_id]
      start_index <- run_start[run_index]
      last_positive_index <- run_end[run_index]
      extinction_observed <- last_positive_index < nrow(dat)
      end_time <- if (extinction_observed) {
        dat$time[last_positive_index + 1L]
      } else {
        NA_real_
      }
      observed_duration <- if (extinction_observed) {
        end_time - dat$time[start_index]
      } else {
        censoring_time - dat$time[start_index]
      }
      established <- observed_duration > threshold
      data.frame(
        combo_id = combo_id,
        replicate = replicate,
        patch = dat$patch[start_index],
        episode_id = episode_id,
        start_time = dat$time[start_index],
        end_time = end_time,
        censoring_time = censoring_time,
        observed_duration = observed_duration,
        extinction_observed = extinction_observed,
        T0 = T0,
        threshold_multiplier = threshold_multiplier,
        fadeout_threshold = threshold,
        survived_beyond_threshold = established,
        early_extinction = extinction_observed && !established,
        fadeout = extinction_observed && established,
        right_censored = !extinction_observed
      )
    }))
  })
  out <- do.call(rbind, episodes)
  if (is.null(out)) {
    out <- data.frame(
      combo_id = character(), replicate = integer(), patch = integer(),
      episode_id = integer(), start_time = numeric(), end_time = numeric(),
      censoring_time = numeric(), observed_duration = numeric(),
      extinction_observed = logical(), T0 = numeric(),
      threshold_multiplier = numeric(), fadeout_threshold = numeric(),
      survived_beyond_threshold = logical(), early_extinction = logical(),
      fadeout = logical(), right_censored = logical()
    )
  }
  rownames(out) <- NULL
  out
}

calculate_occupancy <- function(patch_state, n_patch) {
  patch_state |>
    dplyr::group_by(.data$time) |>
    dplyr::summarise(
      occupied_patches = sum(.data$I > 0),
      occupancy = .data$occupied_patches / n_patch,
      total_I = sum(.data$I),
      .groups = "drop"
    )
}

assign_simulation_year <- function(time, season_period, simulation_end) {
  year <- floor(time / season_period) + 1L
  endpoint_on_boundary <- simulation_end > 0 &&
    abs(simulation_end / season_period -
          round(simulation_end / season_period)) < 1e-9
  if (endpoint_on_boundary) {
    year[abs(time - simulation_end) < 1e-9] <-
      max(1L, as.integer(round(simulation_end / season_period)))
  }
  as.integer(year)
}

detect_global_extinction <- function(occupancy) {
  zero <- occupancy$total_I == 0
  if (!any(zero)) return(NA_real_)
  candidates <- which(zero)
  persistent_zero <- candidates[vapply(
    candidates,
    function(i) all(zero[i:length(zero)]),
    logical(1)
  )]
  if (!length(persistent_zero)) NA_real_ else occupancy$time[persistent_zero[1]]
}

summarise_annual_occupancy <- function(
    occupancy, combo_id, replicate, season_period, dt) {
  simulation_end <- max(occupancy$time)
  expected_observations <- ceiling(season_period / dt)
  occupancy |>
    dplyr::mutate(
      simulation_year = assign_simulation_year(
        .data$time, season_period, simulation_end
      )
    ) |>
    dplyr::group_by(.data$simulation_year) |>
    dplyr::summarise(
      combo_id = combo_id,
      replicate = replicate,
      mean_occupancy = mean(.data$occupancy),
      occupancy_amplitude = max(.data$occupancy) - min(.data$occupancy),
      observations = dplyr::n(),
      complete_year = dplyr::n() >= expected_observations,
      .groups = "drop"
    ) |>
    dplyr::select(
      "combo_id", "replicate", "simulation_year",
      "mean_occupancy", "occupancy_amplitude",
      "observations", "complete_year"
    )
}

fit_occupancy_trend <- function(annual, global_extinction_time, season_period) {
  fit_data <- annual[annual$complete_year, , drop = FALSE]
  fit_through_year <- NA_integer_
  if (is.finite(global_extinction_time)) {
    fit_through_year <- floor(global_extinction_time / season_period) + 1L
    fit_data <- fit_data[
      fit_data$simulation_year <= fit_through_year, , drop = FALSE
    ]
  }
  if (nrow(fit_data) < 2L) {
    return(list(
      intercept = NA_real_, slope = NA_real_, r_squared = NA_real_,
      n_complete_years = nrow(fit_data), fit_through_year = fit_through_year,
      fitted = data.frame()
    ))
  }
  model <- stats::lm(mean_occupancy ~ simulation_year, data = fit_data)
  list(
    intercept = unname(stats::coef(model)[1]),
    slope = unname(stats::coef(model)[2]),
    r_squared = summary(model)$r.squared,
    n_complete_years = nrow(fit_data),
    fit_through_year = fit_through_year,
    fitted = data.frame(
      simulation_year = fit_data$simulation_year,
      fitted_occupancy = stats::predict(model, newdata = fit_data)
    )
  )
}

summarise_fadeout_years <- function(
    episodes, annual, season_period, simulation_end) {
  years <- unique(annual[c("combo_id", "replicate", "simulation_year")])
  fadeouts <- episodes[episodes$fadeout, , drop = FALSE]
  if (nrow(fadeouts)) {
    fadeouts$simulation_year <- assign_simulation_year(
      fadeouts$end_time, season_period, simulation_end
    )
    counts <- fadeouts |>
      dplyr::count(
        .data$combo_id, .data$replicate, .data$simulation_year,
        name = "fadeout_count"
      )
    years <- dplyr::left_join(
      years, counts,
      by = c("combo_id", "replicate", "simulation_year")
    )
  } else {
    years$fadeout_count <- 0L
  }
  years$fadeout_count[is.na(years$fadeout_count)] <- 0L
  years[order(years$combo_id, years$replicate, years$simulation_year), ]
}

make_survival_data <- function(episodes) {
  full <- episodes |>
    dplyr::transmute(
      .data$combo_id, .data$replicate, .data$patch, .data$episode_id,
      survival_time = .data$observed_duration,
      event = as.integer(.data$extinction_observed),
      .data$fadeout_threshold
    )
  established <- episodes |>
    dplyr::filter(.data$survived_beyond_threshold) |>
    dplyr::transmute(
      .data$combo_id, .data$replicate, .data$patch, .data$episode_id,
      survival_time = .data$observed_duration - .data$fadeout_threshold,
      event = as.integer(.data$fadeout)
    )
  list(full = full, established = established)
}

tidy_survfit <- function(fit) {
  summary_fit <- summary(fit)
  if (!length(summary_fit$time)) return(data.frame())
  out <- data.frame(
    time = summary_fit$time,
    survival = summary_fit$surv,
    lower = summary_fit$lower,
    upper = summary_fit$upper
  )
  rbind(
    data.frame(time = 0, survival = 1, lower = 1, upper = 1),
    out
  )
}

plot_annual_diagnostics <- function(annual, trends, parameters, figure_dir) {
  parameters$parameter_label <- vapply(
    seq_len(nrow(parameters)),
    function(i) make_parameter_label(parameters[i, , drop = FALSE]),
    character(1)
  )
  labels <- parameters |>
    dplyr::select("combo_id", "replicate", "parameter_label")
  annual_plot <- annual |>
    dplyr::left_join(
      labels, by = c("combo_id", "replicate")
    ) |>
    dplyr::mutate(
      panel = paste0(.data$parameter_label, " | replicate=", .data$replicate)
    )
  fitted <- trends |>
    dplyr::filter(is.finite(.data$fitted_occupancy)) |>
    dplyr::left_join(labels, by = c("combo_id", "replicate")) |>
    dplyr::mutate(
      panel = paste0(.data$parameter_label, " | replicate=", .data$replicate)
    )

  mean_plot <- ggplot2::ggplot(
    annual_plot,
    ggplot2::aes(.data$simulation_year, .data$mean_occupancy)
  ) +
    ggplot2::geom_line(linewidth = 0.45) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_line(
      data = fitted,
      ggplot2::aes(y = .data$fitted_occupancy),
      colour = "#D55E00", linewidth = 0.8
    ) +
    ggplot2::facet_wrap(~panel, scales = "free_x") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Simulation year", y = "Annual mean occupancy",
      title = "Annual mean patch occupancy",
      caption = "Orange line: exploratory linear fit using complete years only; after global extinction, fit stops in the extinction year."
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 7),
      plot.caption = ggplot2::element_text(hjust = 0)
    )

  amplitude_plot <- ggplot2::ggplot(
    annual_plot,
    ggplot2::aes(.data$simulation_year, .data$occupancy_amplitude)
  ) +
    ggplot2::geom_line(linewidth = 0.45) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::facet_wrap(~panel, scales = "free_x") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Simulation year",
      y = "Annual max occupancy - annual min occupancy",
      title = "Annual occupancy amplitude"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 7))

  ggplot2::ggsave(
    file.path(figure_dir, "annual_mean_occupancy.pdf"),
    mean_plot, width = 11, height = max(5, 3 * ceiling(nrow(labels) / 2))
  )
  ggplot2::ggsave(
    file.path(figure_dir, "annual_occupancy_amplitude.pdf"),
    amplitude_plot, width = 11, height = max(5, 3 * ceiling(nrow(labels) / 2))
  )
}

plot_survival_diagnostics <- function(episodes, figure_dir) {
  survival_data <- make_survival_data(episodes)
  make_plot <- function(dat, conditional = FALSE) {
    if (!nrow(dat)) {
      return(
        ggplot2::ggplot() +
          ggplot2::annotate(
            "text", x = 0, y = 0,
            label = "No eligible episodes"
          ) +
          ggplot2::labs(
            title = if (conditional) {
              "Conditional survival of established infection episodes"
            } else {
              "Survival of all infection episodes"
            }
          ) +
          ggplot2::theme_void() +
          ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0),
            plot.margin = ggplot2::margin(12, 12, 12, 12)
          )
      )
    }
    grouped <- split(dat, interaction(
      dat$combo_id, dat$replicate, drop = TRUE
    ))
    curve <- dplyr::bind_rows(lapply(grouped, function(group) {
      fit <- survival::survfit(
        survival::Surv(survival_time, event) ~ 1, data = group
      )
      cbind(
        data.frame(
          combo_id = group$combo_id[1],
          replicate = group$replicate[1]
        ),
        tidy_survfit(fit)
      )
    }))
    plot <- ggplot2::ggplot(
      curve,
      ggplot2::aes(
        .data$time, .data$survival,
        colour = factor(.data$replicate),
        group = .data$replicate
      )
    ) +
      ggplot2::geom_step(linewidth = 0.6) +
      ggplot2::facet_wrap(~combo_id, scales = "free_x") +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(
        x = if (conditional) {
          "Time since crossing fade-out threshold (days)"
        } else {
          "Time since episode start (days)"
        },
        y = "Episode survival probability",
        colour = "Replicate",
        title = if (conditional) {
          "Conditional survival of established infection episodes"
        } else {
          "Survival of all infection episodes"
        }
      ) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(legend.position = "right")
    if (!conditional) {
      thresholds <- unique(
        episodes[c("combo_id", "fadeout_threshold")]
      )
      plot <- plot + ggplot2::geom_vline(
        data = thresholds,
        ggplot2::aes(xintercept = .data$fadeout_threshold),
        linetype = "dashed", colour = "grey35"
      )
    }
    plot
  }

  ggplot2::ggsave(
    file.path(figure_dir, "episode_survival.pdf"),
    make_plot(survival_data$full), width = 10, height = 6
  )
  ggplot2::ggsave(
    file.path(figure_dir, "established_episode_survival.pdf"),
    make_plot(survival_data$established, conditional = TRUE),
    width = 10, height = 6
  )
}

run_seasonal_fadeout_grid <- function(grid, output_dir, model_file) {
  require_packages(c(
    "dplyr", "ggplot2", "odin", "plagueMetapop", "survival", "tidyr"
  ))
  grid <- validate_seasonal_grid(grid)
  if (!file.exists(model_file)) stop("Seasonal odin model not found: ", model_file)

  data_dir <- file.path(output_dir, "data")
  figure_dir <- file.path(output_dir, "figures")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  generator <- suppressMessages(odin::odin(model_file))

  episode_results <- list()
  annual_results <- list()
  fadeout_year_results <- list()
  summary_results <- list()
  parameter_results <- list()
  extinction_results <- list()
  trend_curve_results <- list()
  index <- 0L

  for (combo_index in seq_len(nrow(grid))) {
    parameters <- grid[combo_index, , drop = FALSE]
    period <- calculate_intrinsic_period(
      parameters$R0, parameters$gamma, parameters$r
    )
    threshold <- parameters$threshold_multiplier * period$T0
    beta <- parameters$R0 * parameters$gamma
    S_star <- parameters$gamma * parameters$K / beta
    I_star <- parameters$r * (parameters$K - S_star) / beta

    for (replicate in seq_len(parameters$n_reps)) {
      index <- index + 1L
      seed <- parameters$base_seed + (combo_index - 1L) * 100000L +
        replicate - 1L
      set.seed(seed)
      S_ini <- stats::rpois(parameters$n_patch, S_star)
      I_ini <- cbind(
        stats::rpois(parameters$n_patch, I_star),
        rep(0L, parameters$n_patch)
      )
      model <- generator$new(
        beta = c(beta, 0),
        gamma = rep(parameters$gamma, 2),
        dt = parameters$dt,
        I_ini = I_ini,
        S_ini = S_ini,
        I2_ini = rep(0L, parameters$n_patch),
        alpha = parameters$alpha,
        strain2_delay = .Machine$integer.max,
        r = rep(parameters$r, parameters$n_patch),
        K = rep(parameters$K, parameters$n_patch),
        season_period = parameters$season_period,
        seasonal_amp = parameters$seasonal_amp,
        peak_day = parameters$peak_day,
        n_patch = parameters$n_patch
      )
      nt <- round(parameters$t_max / parameters$dt)
      raw <- model$run(seq.int(0L, nt))
      if (parameters$dt != 1) raw[, "step"] <- raw[, "step"] * parameters$dt
      runs <- plagueMetapop::conv_odin(raw)
      patch_state <- extract_patch_state(runs)
      plot_occupancy_raster(
        patch_state, parameters, replicate, seed, period$T0, threshold,
        figure_dir
      )
      occupancy <- calculate_occupancy(patch_state, parameters$n_patch)
      episodes <- extract_infection_episodes(
        patch_state, parameters$combo_id, replicate, period$T0,
        parameters$threshold_multiplier
      )
      annual <- summarise_annual_occupancy(
        occupancy, parameters$combo_id, replicate,
        parameters$season_period, parameters$dt
      )
      extinction_time <- detect_global_extinction(occupancy)
      trend <- fit_occupancy_trend(
        annual, extinction_time, parameters$season_period
      )
      annual_fadeouts <- summarise_fadeout_years(
        episodes, annual, parameters$season_period, max(occupancy$time)
      )

      established <- sum(episodes$survived_beyond_threshold)
      fadeouts <- sum(episodes$fadeout)
      patches_with_fadeout <- length(unique(episodes$patch[episodes$fadeout]))
      summary_row <- cbind(
        parameters[c(
          "combo_id", "R0", "gamma", "r", "K", "alpha", "seasonal_amp",
          "season_period", "peak_day", "n_patch", "dt", "t_max",
          "threshold_multiplier"
        )],
        data.frame(
          replicate = replicate,
          seed = seed,
          T0 = period$T0,
          fadeout_threshold = threshold,
          global_extinction = is.finite(extinction_time),
          global_extinction_time = extinction_time,
          total_episodes = nrow(episodes),
          early_extinctions = sum(episodes$early_extinction),
          established_episodes = established,
          observed_fadeouts = fadeouts,
          censored_established_episodes = sum(
            episodes$right_censored &
              episodes$survived_beyond_threshold
          ),
          fadeout_fraction_among_established = if (established > 0) {
            fadeouts / established
          } else {
            NA_real_
          },
          fraction_patches_with_fadeout = patches_with_fadeout /
            parameters$n_patch,
          occupancy_trend_intercept = trend$intercept,
          occupancy_trend_slope = trend$slope,
          occupancy_trend_r_squared = trend$r_squared,
          complete_years_used = trend$n_complete_years,
          trend_fit_through_year = trend$fit_through_year
        )
      )
      parameter_results[[index]] <- cbind(
        parameters, data.frame(replicate = replicate, seed = seed),
        data.frame(
          beta = beta, S_star = S_star, I_star = I_star,
          T0 = period$T0, omega0 = period$omega0,
          eigen_omega0 = period$eigen_omega0,
          fadeout_threshold = threshold
        )
      )
      extinction_results[[index]] <- data.frame(
        combo_id = parameters$combo_id,
        replicate = replicate,
        seed = seed,
        global_extinction = is.finite(extinction_time),
        global_extinction_time = extinction_time
      )
      if (nrow(trend$fitted)) {
        trend_curve_results[[index]] <- cbind(
          data.frame(
            combo_id = parameters$combo_id,
            replicate = replicate
          ),
          trend$fitted
        )
      }
      episode_results[[index]] <- episodes
      annual_results[[index]] <- annual
      fadeout_year_results[[index]] <- annual_fadeouts
      summary_results[[index]] <- summary_row
    }
  }

  episodes <- dplyr::bind_rows(episode_results)
  annual <- dplyr::bind_rows(annual_results)
  annual_fadeouts <- dplyr::bind_rows(fadeout_year_results)
  summaries <- dplyr::bind_rows(summary_results)
  parameters <- dplyr::bind_rows(parameter_results)
  extinctions <- dplyr::bind_rows(extinction_results)
  trends <- dplyr::bind_rows(trend_curve_results)
  if (!ncol(trends)) {
    trends <- data.frame(
      combo_id = character(), replicate = integer(),
      simulation_year = integer(), fitted_occupancy = numeric()
    )
  }

  utils::write.csv(
    episodes, file.path(data_dir, "infection_episodes.csv"), row.names = FALSE
  )
  utils::write.csv(
    annual, file.path(data_dir, "annual_occupancy.csv"), row.names = FALSE
  )
  utils::write.csv(
    annual_fadeouts, file.path(data_dir, "annual_fadeout_counts.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summaries, file.path(data_dir, "replicate_summary.csv"), row.names = FALSE
  )
  utils::write.csv(
    parameters, file.path(data_dir, "parameters_and_seeds.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    extinctions, file.path(data_dir, "global_extinction.csv"), row.names = FALSE
  )

  plot_annual_diagnostics(annual, trends, parameters, figure_dir)
  plot_survival_diagnostics(episodes, figure_dir)
  invisible(list(
    episodes = episodes, annual = annual,
    annual_fadeouts = annual_fadeouts, summaries = summaries,
    parameters = parameters, extinctions = extinctions, trends = trends
  ))
}
