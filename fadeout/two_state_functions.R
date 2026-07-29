## Helpers for the two-state occupancy closure with distinct infected loads in
## persistent and transient source patches.  These are deliberately separate
## from the original one-state approximation.

compute_transient_outbreak_summary <- function(R0, K, r, gamma = 1,
                                                I0 = 10, c = 0.5,
                                                integration_dt = 0.01) {
  values <- c(R0, K, r, gamma, I0, c, integration_dt)
  if (any(!is.finite(values)) || R0 <= 1 || K <= I0 || r <= 0 ||
      gamma <= 0 || I0 <= 0 || c <= 0 || integration_dt <= 0) {
    stop("Invalid deterministic transient-outbreak parameters")
  }

  ## This period closure is provisional.  The implemented population model has
  ## logistic susceptible recruitment rather than constant vital turnover, so
  ## mu=r and rho=gamma are an explicit approximation (as in the earlier
  ## two-state analysis), not an exact parameter identity.
  mu <- r
  rho <- gamma
  exact_argument <- mu * rho * (R0 - 1) - mu^2 * R0^2 / 4
  approximate_argument <- mu * rho * (R0 - 1)
  if (exact_argument <= 0) {
    return(list(
      oscillatory = FALSE, exact_argument = exact_argument,
      T_osc = NA_real_, T_osc_approx = 2 * pi / sqrt(approximate_argument),
      T_T = NA_real_, integral_I = NA_real_, Ibar_T = NA_real_,
      I_peak = NA_real_, time_of_I_peak = NA_real_, trajectory = NULL,
      mapping = "mu approximated by logistic growth parameter r; rho = gamma"
    ))
  }
  T_osc <- 2 * pi / sqrt(exact_argument)
  T_T <- c * T_osc
  times <- sort(unique(c(seq(0, T_T, by = integration_dt), T_T)))
  beta <- R0 * gamma

  ## These are the continuous-time equations corresponding to the repository's
  ## deterministic logistic single-patch model and plagueMetapop::ode_eq():
  ## dS/dt = r*S*(1-S/K) - beta*S*I/K; dI/dt = beta*S*I/K-gamma*I.
  rhs <- function(time, state, parameters) {
    S <- state[["S"]]
    I <- state[["I"]]
    incidence <- beta * S * I / K
    list(c(S = r * S * (1 - S / K) - incidence,
           I = incidence - gamma * I))
  }
  trajectory <- as.data.frame(deSolve::ode(
    y = c(S = K - I0, I = I0), times = times, func = rhs,
    parms = NULL, method = "lsoda", rtol = 1e-10, atol = 1e-12
  ))
  integral_I <- sum(diff(trajectory$time) *
    (head(trajectory$I, -1) + tail(trajectory$I, -1)) / 2)
  Ibar_T <- integral_I / T_T
  peak_index <- which.max(trajectory$I)
  I_peak <- trajectory$I[peak_index]
  time_of_I_peak <- trajectory$time[peak_index]
  if (!is.finite(Ibar_T) || Ibar_T <= 0 || Ibar_T > I_peak + 1e-7) {
    stop("Invalid deterministic transient infected-load summary")
  }
  list(
    oscillatory = TRUE, exact_argument = exact_argument,
    T_osc = T_osc, T_osc_approx = 2 * pi / sqrt(approximate_argument),
    T_T = T_T, integral_I = integral_I, Ibar_T = Ibar_T,
    I_peak = I_peak, time_of_I_peak = time_of_I_peak,
    trajectory = trajectory,
    mapping = "mu approximated by logistic growth parameter r; rho = gamma"
  )
}

compute_transient_outbreak_summary_I1 <- function(
    R0, K, r, I0 = 10, integration_dt = 0.01,
    initial_tmax = 50, maximum_tmax = 800) {
  values <- c(R0, K, r, I0, integration_dt, initial_tmax, maximum_tmax)
  if (any(!is.finite(values)) || R0 <= 1 || K <= I0 || r <= 0 ||
      I0 <= 0 || integration_dt <= 0 || initial_tmax <= 0 ||
      maximum_tmax < initial_tmax) {
    stop("Invalid deterministic hybrid transient-duration parameters")
  }

  ## Disease-generation time normalization: removal rate is exactly 1.
  ## These are the repository's logistic susceptible-growth equations, not a
  ## constant/linear demographic-turnover approximation.
  rhs <- function(time, state, parameters) {
    S <- state[["S"]]
    I <- state[["I"]]
    incidence <- R0 * S * I / K
    list(c(
      S = r * S * (1 - S / K) - incidence,
      I = incidence - I
    ))
  }

  solve_to <- function(tmax) {
    times <- sort(unique(c(seq(0, tmax, by = integration_dt), tmax)))
    as.data.frame(deSolve::ode(
      y = c(S = K - I0, I = I0), times = times, func = rhs,
      parms = NULL, method = "lsoda", rtol = 1e-10, atol = 1e-12
    ))
  }

  interpolate_event <- function(trajectory, index, variable, target) {
    left <- trajectory[index, ]
    right <- trajectory[index + 1L, ]
    left_value <- left[[variable]]
    right_value <- right[[variable]]
    fraction <- (left_value - target) / (left_value - right_value)
    list(
      time = left$time + fraction * (right$time - left$time),
      S = left$S + fraction * (right$S - left$S),
      I = left$I + fraction * (right$I - left$I),
      fraction = fraction,
      left_index = index
    )
  }

  find_events <- function(trajectory) {
    equilibrium_S <- K / R0
    distance <- trajectory$S - equilibrium_S
    interval <- seq_len(nrow(trajectory) - 1L)

    ## dI/dt = I*(R0*S/K - 1). Thus the first downward S=K/R0
    ## crossing is the outbreak peak and the next upward crossing is the first
    ## post-peak trough.
    peak_candidates <- which(
      head(distance, -1L) > 0 & tail(distance, -1L) <= 0
    )
    peak_interval <- if (length(peak_candidates)) {
      peak_candidates[1L]
    } else {
      NA_integer_
    }
    trough_interval <- NA_integer_
    if (!is.na(peak_interval)) {
      trough_candidates <- which(
        interval > peak_interval &
          head(distance, -1L) < 0 & tail(distance, -1L) >= 0
      )
      if (length(trough_candidates)) {
        trough_interval <- trough_candidates[1L]
      }
    }

    peak <- if (!is.na(peak_interval)) {
      interpolate_event(trajectory, peak_interval, "S", equilibrium_S)
    } else {
      NULL
    }
    trough <- if (!is.na(trough_interval)) {
      interpolate_event(trajectory, trough_interval, "S", equilibrium_S)
    } else {
      NULL
    }

    crossing_interval <- NA_integer_
    crossing <- NULL
    if (!is.null(peak) && !is.null(trough)) {
      candidate_crossings <- which(
        head(trajectory$I, -1L) > 1 & tail(trajectory$I, -1L) <= 1 &
          head(trajectory$time, -1L) >= peak$time - integration_dt &
          tail(trajectory$time, -1L) <= trough$time + integration_dt
      )
      if (length(candidate_crossings)) {
        crossing_interval <- candidate_crossings[1L]
        crossing <- interpolate_event(
          trajectory, crossing_interval, "I", 1
        )
      }
    }

    list(
      equilibrium_S = equilibrium_S,
      peak_interval = peak_interval,
      trough_interval = trough_interval,
      crossing_interval = crossing_interval,
      peak = peak,
      trough = trough,
      crossing = crossing
    )
  }

  horizon <- initial_tmax
  trajectory <- NULL
  events <- NULL
  repeat {
    trajectory <- solve_to(horizon)
    events <- find_events(trajectory)
    if (!is.null(events$trough)) break
    if (horizon >= maximum_tmax) break
    horizon <- min(maximum_tmax, horizon * 2)
  }

  peak_found <- !is.null(events$peak)
  trough_found <- !is.null(events$trough)
  crossing_found <- !is.null(events$crossing)

  I_peak <- if (peak_found) events$peak$I else NA_real_
  time_of_I_peak <- if (peak_found) events$peak$time else NA_real_
  S_at_peak <- if (peak_found) events$peak$S else NA_real_
  time_of_first_trough <- if (trough_found) {
    events$trough$time
  } else {
    NA_real_
  }
  I_first_trough <- if (trough_found) events$trough$I else NA_real_
  S_at_first_trough <- if (trough_found) events$trough$S else NA_real_

  endpoint_found <- trough_found
  transient_endpoint <- NA_character_
  T_T <- integral_I <- Ibar_T <- I_at_T_T <- NA_real_
  I_before_crossing <- I_after_crossing <- NA_real_
  crossing_left_time <- crossing_right_time <- NA_real_
  endpoint_left_index <- NA_integer_
  endpoint_tolerance <- 1e-7
  if (endpoint_found) {
    if (I_first_trough <= 1 + endpoint_tolerance) {
      if (crossing_found) {
        transient_endpoint <- "I1_crossing"
        T_T <- events$crossing$time
        I_at_T_T <- events$crossing$I
        endpoint_left_index <- events$crossing$left_index
        I_before_crossing <- trajectory$I[endpoint_left_index]
        I_after_crossing <- trajectory$I[endpoint_left_index + 1L]
        crossing_left_time <- trajectory$time[endpoint_left_index]
        crossing_right_time <- trajectory$time[endpoint_left_index + 1L]
      } else if (abs(I_first_trough - 1) <= endpoint_tolerance) {
        ## Numerical equality at the trough is classified as the I=1 endpoint.
        transient_endpoint <- "I1_crossing"
        T_T <- time_of_first_trough
        I_at_T_T <- 1
        endpoint_left_index <- events$trough$left_index
      } else {
        endpoint_found <- FALSE
      }
    } else {
      transient_endpoint <- "first_trough"
      T_T <- time_of_first_trough
      I_at_T_T <- I_first_trough
      endpoint_left_index <- events$trough$left_index
    }
  }

  if (endpoint_found) {
    before_endpoint <- which(trajectory$time < T_T)
    endpoint_time <- c(trajectory$time[before_endpoint], T_T)
    endpoint_I <- c(trajectory$I[before_endpoint], I_at_T_T)
    if (length(endpoint_time) < 2L) {
      stop("Hybrid transient endpoint does not follow time zero")
    }
    integral_I <- sum(
      diff(endpoint_time) *
        (head(endpoint_I, -1L) + tail(endpoint_I, -1L)) / 2
    )
    Ibar_T <- integral_I / T_T
  }

  status <- if (!trough_found) {
    "no_postpeak_trough"
  } else if (!endpoint_found) {
    "I1_crossing_expected_but_not_found"
  } else if (transient_endpoint == "I1_crossing") {
    "success_I1_crossing"
  } else {
    "success_first_trough"
  }

  ## Old provisional closure, returned only for comparison. It is not exact
  ## for the logistic model and does not define the new T_T.
  old_argument <- r * (R0 - 1) - r^2 * R0^2 / 4
  T_osc_old <- if (old_argument > 0) {
    2 * pi / sqrt(old_argument)
  } else {
    NA_real_
  }

  list(
    endpoint_found = endpoint_found,
    crossing_found = crossing_found,
    status = status,
    T_T = T_T,
    transient_endpoint = transient_endpoint,
    integral_I = integral_I,
    Ibar_T = Ibar_T,
    I_peak = I_peak,
    time_of_I_peak = time_of_I_peak,
    time_of_first_trough = time_of_first_trough,
    I_first_trough = I_first_trough,
    I_at_T_T = I_at_T_T,
    trajectory = trajectory,
    R0 = R0,
    K = K,
    r = r,
    I0 = I0,
    transient_method = "I1_or_first_trough",
    T_osc_old = T_osc_old,
    T_T_old = 0.5 * T_osc_old,
    integration_dt = integration_dt,
    integration_horizon = horizon,
    peak_found = peak_found,
    trough_found = trough_found,
    equilibrium_S = events$equilibrium_S,
    S_at_peak = S_at_peak,
    S_at_first_trough = S_at_first_trough,
    peak_S_left = if (peak_found) {
      trajectory$S[events$peak$left_index]
    } else {
      NA_real_
    },
    peak_S_right = if (peak_found) {
      trajectory$S[events$peak$left_index + 1L]
    } else {
      NA_real_
    },
    trough_S_left = if (trough_found) {
      trajectory$S[events$trough$left_index]
    } else {
      NA_real_
    },
    trough_S_right = if (trough_found) {
      trajectory$S[events$trough$left_index + 1L]
    } else {
      NA_real_
    },
    I_before_crossing = I_before_crossing,
    I_after_crossing = I_after_crossing,
    crossing_left_time = crossing_left_time,
    crossing_right_time = crossing_right_time,
    time_units = "disease_generations"
  )
}

solve_two_state_transient_Ibar <- function(time, p0, q0, P1, alpha,
                                           I_star, Ibar_T, T_T,
                                           numerical_tolerance = 1e-7) {
  inputs <- c(p0, q0, P1, alpha, I_star, Ibar_T, T_T)
  if (any(!is.finite(inputs)) || p0 < 0 || q0 < 0 || p0 + q0 > 1 ||
      P1 < 0 || P1 > 1 || alpha < 0 || I_star <= 0 || Ibar_T <= 0 ||
      T_T <= 0 || any(diff(time) < 0)) {
    stop("Invalid Ibar two-state parameters or initial conditions")
  }
  rhs <- function(t, state, parameters) {
    p <- state[["p"]]
    q <- state[["q"]]
    available <- max(0, 1 - p - q)
    B <- alpha * (I_star * p + Ibar_T * q) * available
    list(c(p = P1 * B, q = (1 - P1) * B - q / T_T))
  }
  ans <- as.data.frame(deSolve::ode(
    y = c(p = p0, q = q0), times = time, func = rhs, parms = NULL,
    method = "lsoda", rtol = 1e-9, atol = 1e-11
  ))
  ans$total <- ans$p + ans$q
  violation <- max(c(-ans$p, -ans$q, ans$total - 1), na.rm = TRUE)
  if (violation > numerical_tolerance) {
    stop("Substantive two-state state-space violation: ", violation)
  }
  ans$p[ans$p < 0] <- 0
  ans$q[ans$q < 0] <- 0
  ans$total <- ans$p + ans$q
  ans
}
