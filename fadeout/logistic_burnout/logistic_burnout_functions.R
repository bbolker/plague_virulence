## Semi-analytical post-epidemic burnout approximation for the normalized
## logistic S-I model. Time is measured in disease generations.

.check_scalar <- function(x, name, lower = -Inf, strict = FALSE) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      if (strict) x <= lower else x < lower) {
    stop(name, " must be one finite numeric value ",
         if (strict) paste0("> ", lower) else paste0(">= ", lower))
  }
}

## Return the endemic equilibrium in densities and host counts.
logistic_endemic_equilibrium <- function(R0, r, K) {
  .check_scalar(R0, "R0", 1, strict = TRUE)
  .check_scalar(r, "r", 0, strict = TRUE)
  .check_scalar(K, "K", 0, strict = TRUE)
  x_star <- 1 / R0
  y_star <- r * (R0 - 1) / R0^2
  list(
    x_star = x_star, y_star = y_star,
    S_star = K * x_star, I_star = K * y_star
  )
}

## Solve the deterministic logistic S-I count model on a fixed output grid.
solve_logistic_SI <- function(R0, r, K, I0 = 1, tmax = 50, dt = 0.02,
                              rtol = 1e-9, atol = 1e-11) {
  .check_scalar(R0, "R0", 1, strict = TRUE)
  .check_scalar(r, "r", 0, strict = TRUE)
  .check_scalar(K, "K", 0, strict = TRUE)
  .check_scalar(I0, "I0", 1)
  .check_scalar(tmax, "tmax", 0, strict = TRUE)
  .check_scalar(dt, "dt", 0, strict = TRUE)
  .check_scalar(rtol, "rtol", 0, strict = TRUE)
  .check_scalar(atol, "atol", 0, strict = TRUE)
  if (I0 >= K) stop("I0 must be smaller than K")
  if (!requireNamespace("deSolve", quietly = TRUE)) {
    stop("Package 'deSolve' is required")
  }

  rhs <- function(time, state, parameters) {
    S <- state[["S"]]
    I <- state[["I"]]
    incidence <- R0 * S * I / K
    list(c(
      S = r * S * (1 - S / K) - incidence,
      I = incidence - I
    ))
  }
  times <- sort(unique(c(seq(0, tmax, by = dt), tmax)))
  out <- as.data.frame(deSolve::ode(
    y = c(S = K - I0, I = I0), times = times,
    func = rhs, parms = NULL, method = "lsoda",
    rtol = rtol, atol = atol
  ))
  if (any(!is.finite(as.matrix(out[c("S", "I")]))) ||
      any(out$S < -sqrt(.Machine$double.eps)) ||
      any(out$I < -sqrt(.Machine$double.eps))) {
    stop("Deterministic logistic S-I solution is non-finite or negative")
  }
  out$S <- pmax(out$S, 0)
  out$I <- pmax(out$I, 0)
  out$x <- out$S / K
  out$y <- out$I / K
  attr(out, "parameters") <- list(
    R0 = R0, r = r, K = K, I0 = I0, dt = dt, tmax = tmax,
    rtol = rtol, atol = atol
  )
  out
}

.interpolate_crossing <- function(trajectory, index, variable, target) {
  left <- trajectory[index, , drop = FALSE]
  right <- trajectory[index + 1L, , drop = FALSE]
  z0 <- left[[variable]]
  z1 <- right[[variable]]
  fraction <- (target - z0) / (z1 - z0)
  list(
    time = left$time + fraction * (right$time - left$time),
    S = left$S + fraction * (right$S - left$S),
    I = left$I + fraction * (right$I - left$I),
    x = left$x + fraction * (right$x - left$x),
    y = left$y + fraction * (right$y - left$y),
    left_index = index,
    fraction = fraction
  )
}

## Locate the first epidemic peak and first post-peak downward y_BL crossing.
find_logistic_boundary_entry <- function(trajectory, R0, r, K,
                                         tolerance = 1e-8) {
  .check_scalar(R0, "R0", 1, strict = TRUE)
  .check_scalar(r, "r", 0, strict = TRUE)
  .check_scalar(K, "K", 0, strict = TRUE)
  required <- c("time", "S", "I", "x", "y")
  if (!is.data.frame(trajectory) || nrow(trajectory) < 3L ||
      any(!required %in% names(trajectory))) {
    stop("trajectory must contain time, S, I, x, and y")
  }
  eq <- logistic_endemic_equilibrium(R0, r, K)
  y_BL <- eq$y_star
  I_BL <- K * y_BL
  n <- nrow(trajectory)

  ## Since dI/dt=(R0*x-1)I, the first downward x=1/R0 crossing
  ## is the first epidemic peak.
  peak_intervals <- which(
    head(trajectory$x, -1L) > eq$x_star &
      tail(trajectory$x, -1L) <= eq$x_star
  )
  if (!length(peak_intervals)) {
    return(list(
      entry_found = FALSE, status = "no_epidemic_peak",
      t_peak = NA_real_, I_peak = NA_real_, t_in = NA_real_,
      x_in = NA_real_, y_in = NA_real_, S_in = NA_real_, I_in = NA_real_,
      y_BL = y_BL, I_BL = I_BL, dy_dt_in = NA_real_
    ))
  }
  peak <- .interpolate_crossing(
    trajectory, peak_intervals[1L], "x", eq$x_star
  )
  interval <- seq_len(n - 1L)
  entry_intervals <- which(
    interval >= peak$left_index &
      head(trajectory$y, -1L) > y_BL &
      tail(trajectory$y, -1L) <= y_BL
  )
  if (!length(entry_intervals)) {
    return(list(
      entry_found = FALSE, status = "no_postpeak_downward_boundary_crossing",
      t_peak = peak$time, I_peak = peak$I, t_in = NA_real_,
      x_in = NA_real_, y_in = NA_real_, S_in = NA_real_, I_in = NA_real_,
      y_BL = y_BL, I_BL = I_BL, dy_dt_in = NA_real_
    ))
  }
  entry <- .interpolate_crossing(
    trajectory, entry_intervals[1L], "y", y_BL
  )
  dy_dt_in <- (R0 * entry$x - 1) * entry$y
  valid <- entry$time > peak$time &&
    dy_dt_in < tolerance &&
    abs(entry$y - y_BL) <= tolerance * max(1, y_BL)
  list(
    entry_found = valid,
    status = if (valid) "success" else "invalid_interpolated_entry",
    t_peak = peak$time, I_peak = peak$I,
    t_in = entry$time, x_in = entry$x, y_in = entry$y,
    S_in = entry$S, I_in = entry$I,
    y_BL = y_BL, I_BL = I_BL, dy_dt_in = dy_dt_in,
    peak_left_index = peak$left_index,
    entry_left_index = entry$left_index
  )
}

## Logistic susceptible recovery after boundary-layer entry.
logistic_recovery_x <- function(t, x_in, r) {
  if (!is.numeric(t) || any(!is.finite(t)) || any(t < 0)) {
    stop("t must contain finite non-negative values")
  }
  .check_scalar(x_in, "x_in", 0, strict = TRUE)
  if (x_in > 1) stop("x_in must be <= 1")
  .check_scalar(r, "r", 0, strict = TRUE)
  stats::plogis(stats::qlogis(x_in) + r * t)
}

.logspace_add <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log(exp(a - m) + exp(b - m))
  both_negative_infinite <- is.infinite(a) & a < 0 &
    is.infinite(b) & b < 0
  out[both_negative_infinite] <- -Inf
  out
}

## Integrated per-lineage net growth H(t)=int_0^t [R0*x(s)-1] ds.
logistic_lineage_H <- function(t, R0, r, x_in) {
  if (!is.numeric(t) || any(!is.finite(t)) || any(t < 0)) {
    stop("t must contain finite non-negative values")
  }
  .check_scalar(R0, "R0", 1, strict = TRUE)
  .check_scalar(r, "r", 0, strict = TRUE)
  .check_scalar(x_in, "x_in", 0, strict = TRUE)
  if (x_in > 1) stop("x_in must be <= 1")
  log_A <- if (x_in == 1) -Inf else log1p(-x_in) - log(x_in)
  log_numerator <- .logspace_add(r * t, log_A)
  log_denominator <- .logspace_add(0, log_A)
  (R0 / r) * (log_numerator - log_denominator) - t
}

## Kendall extinction probability for one lineage in logistic recovery.
logistic_lineage_extinction_probability <- function(
    R0, r, x_in, rel.tol = 1e-8, abs.tol = 0,
    subdivisions = 500L) {
  .check_scalar(R0, "R0", 1, strict = TRUE)
  .check_scalar(r, "r", 0, strict = TRUE)
  .check_scalar(x_in, "x_in", 0, strict = TRUE)
  if (x_in > 1) stop("x_in must be <= 1")
  .check_scalar(rel.tol, "rel.tol", 0, strict = TRUE)
  .check_scalar(abs.tol, "abs.tol", 0)
  if (length(subdivisions) != 1L || !is.finite(subdivisions) ||
      subdivisions < 50) {
    stop("subdivisions must be at least 50")
  }

  warnings <- character()
  ans <- tryCatch(
    withCallingHandlers(
      stats::integrate(
        function(t) exp(-logistic_lineage_H(t, R0, r, x_in)),
        lower = 0, upper = Inf, rel.tol = rel.tol, abs.tol = abs.tol,
        subdivisions = as.integer(subdivisions), stop.on.error = FALSE
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(ans, "error")) {
    return(list(
      J = NA_real_, q1 = NA_real_, absolute.error = NA_real_,
      subdivisions = NA_integer_, message = conditionMessage(ans),
      warnings = paste(unique(warnings), collapse = " | "),
      converged = FALSE
    ))
  }
  J <- ans$value
  q1 <- if (is.finite(J) && J >= 0) stats::plogis(log(J)) else NA_real_
  converged <- is.finite(J) && J >= 0 && is.finite(q1) &&
    q1 >= 0 && q1 <= 1 && identical(ans$message, "OK")
  list(
    J = J, q1 = q1, absolute.error = ans$abs.error,
    subdivisions = ans$subdivisions,
    message = ans$message,
    warnings = paste(unique(warnings), collapse = " | "),
    converged = converged
  )
}

## Full deterministic-entry plus branching-process burnout approximation.
logistic_burnout_probability <- function(
    R0, r, K, I0 = 1,
    lineage_count_method = c("round", "continuous", "floor", "ceiling"),
    initial_tmax = 50, maximum_tmax = 1600, dt = 0.02,
    rtol = 1e-9, atol = 1e-11,
    rel.tol = 1e-8, abs.tol = 0, subdivisions = 500L,
    keep_trajectory = FALSE) {
  lineage_count_method <- match.arg(lineage_count_method)
  .check_scalar(R0, "R0", 1, strict = TRUE)
  .check_scalar(r, "r", 0, strict = TRUE)
  .check_scalar(K, "K", 0, strict = TRUE)
  .check_scalar(I0, "I0", 1)
  .check_scalar(initial_tmax, "initial_tmax", 0, strict = TRUE)
  .check_scalar(maximum_tmax, "maximum_tmax", initial_tmax)
  if (I0 >= K) stop("I0 must be smaller than K")

  eq <- logistic_endemic_equilibrium(R0, r, K)
  horizon <- initial_tmax
  entry <- NULL
  trajectory <- NULL
  deterministic_message <- "OK"
  repeat {
    trajectory <- tryCatch(
      solve_logistic_SI(
        R0, r, K, I0, horizon, dt = dt, rtol = rtol, atol = atol
      ),
      error = function(e) e
    )
    if (inherits(trajectory, "error")) {
      deterministic_message <- conditionMessage(trajectory)
      break
    }
    entry <- find_logistic_boundary_entry(trajectory, R0, r, K)
    if (entry$entry_found ||
        entry$status == "no_epidemic_peak" ||
        horizon >= maximum_tmax) break
    horizon <- min(maximum_tmax, 2 * horizon)
  }

  base <- list(
    R0 = R0, r = r, K = K, I0 = I0,
    x_star = eq$x_star, y_star = eq$y_star,
    t_peak = if (is.null(entry)) NA_real_ else entry$t_peak,
    I_peak = if (is.null(entry)) NA_real_ else entry$I_peak,
    t_in = if (is.null(entry)) NA_real_ else entry$t_in,
    x_in = if (is.null(entry)) NA_real_ else entry$x_in,
    y_in = if (is.null(entry)) NA_real_ else entry$y_in,
    I_in = if (is.null(entry)) NA_real_ else entry$I_in,
    y_BL = eq$y_star,
    I_BL = K * eq$y_star,
    m_raw = K * eq$y_star,
    m_used = NA_real_, q1 = NA_real_, J = NA_real_,
    log_P_burnout = NA_real_, P_burnout = NA_real_,
    lineage_count_method = lineage_count_method,
    status = if (is.null(entry)) "deterministic_solver_failure" else entry$status,
    deterministic_converged = !inherits(trajectory, "error"),
    boundary_entry_found = !is.null(entry) && entry$entry_found,
    integration_converged = FALSE,
    deterministic_message = deterministic_message,
    integration_message = NA_character_,
    integration_absolute_error = NA_real_,
    integration_subdivisions = NA_integer_,
    integration_horizon = horizon,
    trajectory = if (keep_trajectory &&
      !inherits(trajectory, "error")) trajectory else NULL
  )
  if (!base$boundary_entry_found) return(base)
  if (!is.finite(entry$x_in) || entry$x_in <= 0 || entry$x_in > 1) {
    base$status <- "invalid_x_in"
    return(base)
  }

  extinction <- logistic_lineage_extinction_probability(
    R0, r, entry$x_in, rel.tol, abs.tol, subdivisions
  )
  base$q1 <- extinction$q1
  base$J <- extinction$J
  base$integration_converged <- extinction$converged
  base$integration_message <- paste(
    c(extinction$message, extinction$warnings[nzchar(extinction$warnings)]),
    collapse = " | "
  )
  base$integration_absolute_error <- extinction$absolute.error
  base$integration_subdivisions <- extinction$subdivisions
  if (!extinction$converged) {
    base$status <- "integration_failure"
    return(base)
  }

  m_raw <- base$m_raw
  m_used <- switch(
    lineage_count_method,
    round = max(1, round(m_raw)),
    floor = max(1, floor(m_raw)),
    ceiling = max(1, ceiling(m_raw)),
    continuous = m_raw
  )
  log_P <- m_used * log(extinction$q1)
  base$m_used <- m_used
  base$log_P_burnout <- log_P
  base$P_burnout <- exp(log_P)
  if (!is.finite(base$P_burnout) ||
      base$P_burnout < 0 || base$P_burnout > 1) {
    base$status <- "invalid_burnout_probability"
    return(base)
  }
  base$status <- if (m_raw < 1) {
    paste0("success_boundary_count_below_one_", lineage_count_method,
           "_minimum_one")
  } else {
    "success"
  }
  base
}

## Exact thinning simulation for the nonhomogeneous linear birth-death process.
simulate_logistic_lineage <- function(
    R0, r, x_in, n_replicates = 1000L, seed = 1L,
    tmax = 300, population_cap = 500L) {
  .check_scalar(R0, "R0", 1, strict = TRUE)
  .check_scalar(r, "r", 0, strict = TRUE)
  .check_scalar(x_in, "x_in", 0, strict = TRUE)
  if (x_in > 1) stop("x_in must be <= 1")
  if (length(n_replicates) != 1L || n_replicates < 1) {
    stop("n_replicates must be positive")
  }
  .check_scalar(tmax, "tmax", 0, strict = TRUE)
  if (length(population_cap) != 1L || population_cap < 2) {
    stop("population_cap must be at least 2")
  }
  set.seed(seed)
  extinct <- logical(n_replicates)
  terminal_I <- integer(n_replicates)
  terminal_time <- numeric(n_replicates)
  envelope <- R0 + 1
  for (replicate in seq_len(n_replicates)) {
    time <- 0
    I <- 1L
    while (I > 0L && I < population_cap && time < tmax) {
      time <- time + stats::rexp(1L, rate = envelope * I)
      if (time >= tmax) break
      lambda <- R0 * logistic_recovery_x(time, x_in, r)
      if (stats::runif(1L) <= (lambda + 1) / envelope) {
        if (stats::runif(1L) <= lambda / (lambda + 1)) {
          I <- I + 1L
        } else {
          I <- I - 1L
        }
      }
    }
    extinct[replicate] <- I == 0L
    terminal_I[replicate] <- I
    terminal_time[replicate] <- time
  }
  list(
    extinction_frequency = mean(extinct),
    monte_carlo_se = sqrt(mean(extinct) * (1 - mean(extinct)) / n_replicates),
    n_replicates = n_replicates,
    n_extinct = sum(extinct),
    n_capped = sum(terminal_I >= population_cap),
    n_time_censored = sum(terminal_I > 0 & terminal_I < population_cap),
    tmax = tmax, population_cap = population_cap,
    seed = seed, terminal_time = terminal_time
  )
}
