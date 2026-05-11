## Two-strain one-wave burnout approximation with ODE boundary detection.
##
## Main function:
##   two_strain_wave_burnout()
##
## The calculation is conditional on using the deterministic ODE trajectory after
## initial establishment. It does not multiply by early fizzle probabilities.
##
## Defaults:
##   y01 = k1 / N
##   y02 = k2 / N
##   x0  = 1 - y01 - y02
##
## Notation follows the derivation as closely as possible:
##   ystar_i  boundary-layer prevalence y_i^*
##   xa       susceptible fraction when the first strain enters the boundary layer
##   xb       susceptible fraction when the second strain enters the boundary layer
##   q_i_C    one-infective extinction probability in the common boundary layer
##   G        rare/resident generating function, e.g. G_{1|2}(q_1^C)
##   Q_i      one-wave burnout probability for strain i

## -------------------------------------------------------------------------
## Single-strain burnout interface
## -------------------------------------------------------------------------

## Boundary-layer prevalence y^* = epsilon * (1 - 1/R0).
ystar <- function(R0, epsilon) {
  if (any(R0 <= 1)) stop("R0 must be > 1.")
  if (epsilon <= 0) stop("epsilon must be positive.")
  epsilon * (1 - 1 / R0)
}

## Common-boundary one-infective extinction probability q_i^C.
## By default this calls burnout::q_approx(), matching the burnout package.
## A custom q_fun may be supplied if desired; it should accept named arguments
## R0, epsilon, and xin.
qC_single <- function(R0, epsilon, xb, q_fun = NULL,
                      q_method = c("q_approx", "q_exact")) {
  if (!is.null(q_fun)) {
    return(q_fun(R0 = R0, epsilon = epsilon, xin = xb))
  }

  if (!requireNamespace("burnout", quietly = TRUE)) {
    stop("Please install the burnout package: remotes::install_github('davidearn/burnout')",
         call. = FALSE)
  }

  q_method <- match.arg(q_method)
  q_fun_burnout <- getExportedValue("burnout", q_method)
  q_fun_burnout(R0 = R0, epsilon = epsilon, xin = xb)
}

## -------------------------------------------------------------------------
## Rare/resident generating function for sequential boundary entry
## -------------------------------------------------------------------------

## Resident trajectory after the first strain enters its boundary layer.
## If strain 1 enters first, this is y_2(x):
##   y_2(x) = y_{2,a} + x_a - x + (1/R02) log(x/x_a).
Y_late <- function(x, R0_late, xa, y_late_a) {
  y_late_a + xa - x + (1 / R0_late) * log(x / xa)
}

## Solve for xb, the susceptible fraction when the late-entering strain reaches
## its own boundary layer.
solve_xb <- function(R0_late, xa, y_late_a, ystar_late) {
  if (y_late_a <= ystar_late) {
    stop("The late-entering strain is already at or below its boundary layer at xa.")
  }
  if (!requireNamespace("gsl", quietly = TRUE)) {
    stop("Please install the gsl package: install.packages('gsl')", call. = FALSE)
  }

  C <- y_late_a + xa - ystar_late
  W_arg <- -R0_late * xa * exp(-R0_late * C)

  if (W_arg < -1 / exp(1) || W_arg >= 0) {
    stop("Lambert-W argument is outside the real principal-branch domain.")
  }

  xb <- -as.numeric(gsl::lambert_W0(W_arg)) / R0_late

  if (!is.finite(xb) || xb <= 0 || xb >= xa) {
    stop("Lambert-W formula returned an invalid xb. Check the entry state.")
  }

  xb
}

## Trapezoid integral on a sorted x_grid.
trapz <- function(x_grid, y_grid) {
  if (length(x_grid) != length(y_grid)) stop("x_grid and y_grid must have the same length.")
  if (length(x_grid) < 2L) stop("Need at least two grid points.")
  sum(diff(x_grid) * (head(y_grid, -1) + tail(y_grid, -1)) / 2)
}

## Cumulative trapezoid integral from x_grid[1] to x_grid[i].
cumtrapz <- function(x_grid, y_grid) {
  out <- numeric(length(x_grid))
  if (length(x_grid) == 1L) return(out)
  increments <- diff(x_grid) * (head(y_grid, -1) + tail(y_grid, -1)) / 2
  out[-1] <- cumsum(increments)
  out
}

## Compute the rare/resident generating function.
## If strain 1 enters first and strain 2 enters later, this is G_{1|2}(z).
## If strain 2 enters first, labels are swapped by the caller.
G_first_given_late <- function(z, R0_first, R0_late, xa, y_late_a, xb,
                               n_x = 4000L) {
  if (z < 0 || z > 1) stop("z must be in [0, 1].")
  if (xb <= 0 || xa <= xb) stop("Need 0 < xb < xa.")
  if (n_x < 50L) warning("n_x is small; accuracy may be poor.")

  x_grid <- seq(xb, xa, length.out = n_x)
  y_late <- Y_late(x_grid, R0_late = R0_late, xa = xa, y_late_a = y_late_a)

  if (any(!is.finite(y_late)) || any(y_late <= 0)) {
    stop("Late-strain prevalence became non-positive on [xb, xa]. Check xb/xa/state.")
  }

  ## h(u) = (R0_first*u - 1) / (R0_late*u*y_late(u)).
  h_grid <- (R0_first * x_grid - 1) / (R0_late * x_grid * y_late)

  ## H(u) = integral_u^xa h(v) dv.
  H_from_xb <- cumtrapz(x_grid, h_grid)
  H_total <- tail(H_from_xb, 1)
  H_grid <- H_total - H_from_xb
  H_xb <- H_grid[1]

  ## A = integral_xb^xa [R0_first / (R0_late*y_late(u))] * exp(H(u)) du.
  A_integrand <- (R0_first / (R0_late * y_late)) * exp(H_grid)
  if (any(!is.finite(A_integrand))) {
    stop("The A integral overflowed. A log-scale implementation may be needed for these parameters.")
  }
  A <- trapz(x_grid, A_integrand)

  if (z == 1) {
    G <- 1
  } else {
    G <- 1 - exp(H_xb) * (1 - z) / (1 + (1 - z) * A)
  }

  G <- min(max(G, 0), 1)

  list(
    G = G,
    A = A,
    H_xb = H_xb,
    x_grid = x_grid,
    H_grid = H_grid,
    y_late_grid = y_late
  )
}

## Sequential-entry approximation: first strain enters its boundary layer,
## late strain remains deterministic and enters later.
sequential_entry_wave <- function(R0_first, R0_late, epsilon, N,
                                  xa, y_late_a,
                                  ystar_first = NULL,
                                  ystar_late = NULL,
                                  q_fun = NULL,
                                  q_method = c("q_approx", "q_exact"),
                                  n_x = 4000L,
                                  discrete_floor = FALSE) {
  q_method <- match.arg(q_method)

  if (is.null(ystar_first)) ystar_first <- ystar(R0_first, epsilon)
  if (is.null(ystar_late)) ystar_late <- ystar(R0_late, epsilon)

  if (discrete_floor) {
    ystar_first <- max(ystar_first, 1 / N)
    ystar_late <- max(ystar_late, 1 / N)
  }

  m_first <- N * ystar_first
  m_late <- N * ystar_late

  xb <- solve_xb(
    R0_late = R0_late,
    xa = xa,
    y_late_a = y_late_a,
    ystar_late = ystar_late
  )

  q_first_C <- qC_single(R0 = R0_first, epsilon = epsilon, xb = xb,
                         q_fun = q_fun, q_method = q_method)
  q_late_C <- qC_single(R0 = R0_late, epsilon = epsilon, xb = xb,
                        q_fun = q_fun, q_method = q_method)

  G_out <- G_first_given_late(
    z = q_first_C,
    R0_first = R0_first,
    R0_late = R0_late,
    xa = xa,
    y_late_a = y_late_a,
    xb = xb,
    n_x = n_x
  )

  Q_first <- G_out$G^m_first
  Q_late <- q_late_C^m_late

  outcomes_first_late <- c(
    both_burnout = Q_first * Q_late,
    first_only_burnout = Q_first * (1 - Q_late),
    late_only_burnout = (1 - Q_first) * Q_late,
    neither_burnout = (1 - Q_first) * (1 - Q_late)
  )

  list(
    branch = "sequential_entry",
    ystar = c(first = ystar_first, late = ystar_late),
    m = c(first = m_first, late = m_late),
    xa = xa,
    xb = xb,
    y_late_a = y_late_a,
    qC = c(first = q_first_C, late = q_late_C),
    G_first_given_late = list(
      G_at_qC = G_out$G,
      A = G_out$A,
      H_xb = G_out$H_xb
    ),
    burnout = c(first = Q_first, late = Q_late),
    outcomes_first_late = outcomes_first_late
  )
}

## Common-boundary approximation: start both strains directly in the common boundary layer.
common_boundary_wave <- function(R01, R02, epsilon, N, xb, m1, m2,
                                 q_fun = NULL,
                                 q_method = c("q_approx", "q_exact")) {
  q_method <- match.arg(q_method)

  m1 <- max(m1, 0)
  m2 <- max(m2, 0)

  q1_C <- qC_single(R0 = R01, epsilon = epsilon, xb = xb,
                    q_fun = q_fun, q_method = q_method)
  q2_C <- qC_single(R0 = R02, epsilon = epsilon, xb = xb,
                    q_fun = q_fun, q_method = q_method)

  Q1 <- q1_C^m1
  Q2 <- q2_C^m2

  outcomes_by_strain <- c(
    both_burnout = Q1 * Q2,
    strain1_only_burnout = Q1 * (1 - Q2),
    strain2_only_burnout = (1 - Q1) * Q2,
    neither_burnout = (1 - Q1) * (1 - Q2)
  )

  list(
    branch = "common_boundary",
    xb = xb,
    m = c(m1 = m1, m2 = m2),
    qC = c(q1_C = unname(q1_C), q2_C = unname(q2_C)),
    qC_by_strain = c(q1_C = unname(q1_C), q2_C = unname(q2_C)),
    strain_burnout = c(Q1 = unname(Q1), Q2 = unname(Q2)),
    outcomes_by_strain = outcomes_by_strain
  )
}

## -------------------------------------------------------------------------
## Deterministic two-strain ODE and boundary-event detection
## -------------------------------------------------------------------------

## ODE state is c(x, y1, y2). Time is in units of the mean infectious period.
two_strain_rhs <- function(state, R01, R02, epsilon) {
  x <- state["x"]
  y1 <- state["y1"]
  y2 <- state["y2"]

  dx <- epsilon * (1 - x) - x * (R01 * y1 + R02 * y2)
  dy1 <- (R01 * x - 1) * y1
  dy2 <- (R02 * x - 1) * y2

  c(x = dx, y1 = dy1, y2 = dy2)
}

rk4_step <- function(state, dt, R01, R02, epsilon) {
  k_1 <- two_strain_rhs(state, R01, R02, epsilon)
  k_2 <- two_strain_rhs(state + 0.5 * dt * k_1, R01, R02, epsilon)
  k_3 <- two_strain_rhs(state + 0.5 * dt * k_2, R01, R02, epsilon)
  k_4 <- two_strain_rhs(state + dt * k_3, R01, R02, epsilon)
  state + dt * (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6
}

linear_event_state <- function(prev, state_next, alpha) {
  prev + alpha * (state_next - prev)
}

## Run ODE until a post-escape downward crossing of a boundary threshold is detected.
run_two_strain_ode_until_boundary <- function(R01, R02, epsilon, N,
                                              k1 = 1, k2 = 1,
                                              y01 = k1 / N,
                                              y02 = k2 / N,
                                              x0 = 1 - y01 - y02,
                                              ystar1 = NULL,
                                              ystar2 = NULL,
                                              t_max = 200,
                                              dt = 0.01,
                                              discrete_floor = FALSE,
                                              store_trajectory = TRUE) {
  if (is.null(ystar1)) ystar1 <- ystar(R01, epsilon)
  if (is.null(ystar2)) ystar2 <- ystar(R02, epsilon)

  if (discrete_floor) {
    ystar1 <- max(ystar1, 1 / N)
    ystar2 <- max(ystar2, 1 / N)
  }

  if (x0 <= 0 || y01 < 0 || y02 < 0) stop("Invalid initial state.")
  if (dt <= 0 || t_max <= 0) stop("dt and t_max must be positive.")

  state <- c(x = x0, y1 = y01, y2 = y02)
  time <- 0

  ## If an initial prevalence is already above its threshold, treat it as escaped.
  escaped <- c(strain1 = y01 >= ystar1, strain2 = y02 >= ystar2)

  n_steps <- ceiling(t_max / dt)

  traj <- NULL
  if (store_trajectory) {
    traj <- matrix(NA_real_, nrow = n_steps + 1L, ncol = 6L)
    colnames(traj) <- c("time", "x", "y1", "y2", "escaped1", "escaped2")
    traj[1, ] <- c(time, state, as.numeric(escaped))
    traj_i <- 1L
  }

  for (step in seq_len(n_steps)) {
    prev <- state
    prev_time <- time
    next_state <- rk4_step(prev, dt, R01, R02, epsilon)
    next_state <- pmax(next_state, 0)
    names(next_state) <- c("x", "y1", "y2")
    time <- prev_time + dt

    ## Upward crossings: leaving the initial low-infection region.
    if (!escaped["strain1"] && prev["y1"] < ystar1 && next_state["y1"] >= ystar1) {
      escaped["strain1"] <- TRUE
    }
    if (!escaped["strain2"] && prev["y2"] < ystar2 && next_state["y2"] >= ystar2) {
      escaped["strain2"] <- TRUE
    }

    ## Downward crossings after escape: entering the burnout boundary layer.
    alpha_down <- c(strain1 = Inf, strain2 = Inf)
    if (escaped["strain1"] && prev["y1"] > ystar1 && next_state["y1"] <= ystar1) {
      alpha_down["strain1"] <- (prev["y1"] - ystar1) / (prev["y1"] - next_state["y1"])
    }
    if (escaped["strain2"] && prev["y2"] > ystar2 && next_state["y2"] <= ystar2) {
      alpha_down["strain2"] <- (prev["y2"] - ystar2) / (prev["y2"] - next_state["y2"])
    }

    if (any(is.finite(alpha_down))) {
      first_alpha <- min(alpha_down)
      event_state <- linear_event_state(prev, next_state, first_alpha)
      event_time <- prev_time + first_alpha * dt
      first_strain <- names(which.min(alpha_down))

      ## If both downward crossings occur in the same RK4 step, treat this as common-boundary entry.
      same_step <- is.finite(alpha_down) & abs(alpha_down - first_alpha) <= 1e-6

      event <- list(
        boundary_event_found = TRUE,
        time = event_time,
        state = event_state,
        first_strain = first_strain,
        escaped = escaped,
        same_step_down = same_step,
        ystar = c(ystar1 = ystar1, ystar2 = ystar2),
        alpha_down = alpha_down
      )

      if (store_trajectory) {
        traj <- traj[seq_len(traj_i), , drop = FALSE]
        event$trajectory <- as.data.frame(traj)
      }
      return(event)
    }

    state <- next_state

    if (store_trajectory) {
      traj_i <- traj_i + 1L
      traj[traj_i, ] <- c(time, state, as.numeric(escaped))
    }
  }

  out <- list(
    boundary_event_found = FALSE,
    time = time,
    state = state,
    escaped = escaped,
    ystar = c(ystar1 = ystar1, ystar2 = ystar2)
  )
  if (store_trajectory) {
    traj <- traj[seq_len(traj_i), , drop = FALSE]
    out$trajectory <- as.data.frame(traj)
  }
  out
}

## -------------------------------------------------------------------------
## Main function
## -------------------------------------------------------------------------

## One-wave two-strain burnout probability calculator.
##
## The result is conditional on using the ODE after initial establishment.
## It does not multiply by probabilities of avoiding early fizzle.
two_strain_wave_burnout <- function(R01, R02, epsilon, N,
                                         k1 = 1, k2 = 1,
                                         y01 = k1 / N,
                                         y02 = k2 / N,
                                         x0 = 1 - y01 - y02,
                                         ystar1 = NULL,
                                         ystar2 = NULL,
                                         q_fun = NULL,
                                         q_method = c("q_approx", "q_exact"),
                                         n_x = 4000L,
                                         discrete_floor = FALSE,
                                         t_max = 200,
                                         dt = 0.01,
                                         store_trajectory = TRUE) {
  q_method <- match.arg(q_method)

  if (is.null(ystar1)) ystar1 <- ystar(R01, epsilon)
  if (is.null(ystar2)) ystar2 <- ystar(R02, epsilon)

  if (discrete_floor) {
    ystar1 <- max(ystar1, 1 / N)
    ystar2 <- max(ystar2, 1 / N)
  }

  ode_event <- run_two_strain_ode_until_boundary(
    R01 = R01,
    R02 = R02,
    epsilon = epsilon,
    N = N,
    k1 = k1,
    k2 = k2,
    y01 = y01,
    y02 = y02,
    x0 = x0,
    ystar1 = ystar1,
    ystar2 = ystar2,
    t_max = t_max,
    dt = dt,
    discrete_floor = FALSE,  ## already applied above if requested
    store_trajectory = store_trajectory
  )

  inputs <- list(
    R01 = R01, R02 = R02, epsilon = epsilon, N = N,
    k1 = k1, k2 = k2, y01 = y01, y02 = y02, x0 = x0,
    t_max = t_max, dt = dt, q_method = q_method
  )

  if (!isTRUE(ode_event$boundary_event_found)) {
    return(list(
      boundary_event_found = FALSE,
      boundary_event = "no",
      message = "No post-escape downward boundary crossing was detected before t_max.",
      inputs = inputs,
      ystar = c(ystar1 = ystar1, ystar2 = ystar2),
      ode = ode_event
    ))
  }

  event_state <- ode_event$state
  x_event <- unname(event_state["x"])
  y1_event <- unname(event_state["y1"])
  y2_event <- unname(event_state["y2"])

  first_strain <- ode_event$first_strain
  escaped <- ode_event$escaped
  same_step <- ode_event$same_step_down

  ## Branch A: both strains re-enter in the same RK4 step -> common boundary layer.
  if (all(same_step)) {
    common <- common_boundary_wave(
      R01 = R01, R02 = R02, epsilon = epsilon, N = N,
      xb = x_event,
      m1 = N * ystar1,
      m2 = N * ystar2,
      q_fun = q_fun,
      q_method = q_method
    )

    return(c(common, list(
      boundary_event_found = TRUE,
      boundary_event = "yes",
      entry_case = "both_strains_reentered_same_step",
      first_strain = first_strain,
      event_time = ode_event$time,
      event_state = event_state,
      ystar = c(ystar1 = ystar1, ystar2 = ystar2),
      inputs = inputs,
      ode = ode_event
    )))
  }

  ## Branch B: one strain re-enters before the other ever escapes.
  ## Treat event time as common-boundary entry. The re-entering strain uses
  ## its threshold count; the non-escaped strain uses its current ODE count.
  if (first_strain == "strain1" && !escaped["strain2"]) {
    common <- common_boundary_wave(
      R01 = R01, R02 = R02, epsilon = epsilon, N = N,
      xb = x_event,
      m1 = N * ystar1,
      m2 = N * y2_event,
      q_fun = q_fun,
      q_method = q_method
    )

    return(c(common, list(
      boundary_event_found = TRUE,
      boundary_event = "yes",
      entry_case = "strain1_reentered_before_strain2_escaped",
      first_strain = "strain1",
      event_time = ode_event$time,
      event_state = event_state,
      ystar = c(ystar1 = ystar1, ystar2 = ystar2),
      inputs = inputs,
      ode = ode_event
    )))
  }

  if (first_strain == "strain2" && !escaped["strain1"]) {
    common <- common_boundary_wave(
      R01 = R01, R02 = R02, epsilon = epsilon, N = N,
      xb = x_event,
      m1 = N * y1_event,
      m2 = N * ystar2,
      q_fun = q_fun,
      q_method = q_method
    )

    return(c(common, list(
      boundary_event_found = TRUE,
      boundary_event = "yes",
      entry_case = "strain2_reentered_before_strain1_escaped",
      first_strain = "strain2",
      event_time = ode_event$time,
      event_state = event_state,
      ystar = c(ystar1 = ystar1, ystar2 = ystar2),
      inputs = inputs,
      ode = ode_event
    )))
  }

  ## Branch C: sequential entry after both strains have escaped.
  if (first_strain == "strain1") {
    ## strain 1 enters first, strain 2 is late/resident.
    seq_out <- sequential_entry_wave(
      R0_first = R01,
      R0_late = R02,
      epsilon = epsilon,
      N = N,
      xa = x_event,
      y_late_a = y2_event,
      ystar_first = ystar1,
      ystar_late = ystar2,
      q_fun = q_fun,
      q_method = q_method,
      n_x = n_x,
      discrete_floor = FALSE
    )

    Q1 <- unname(seq_out$burnout["first"])
    Q2 <- unname(seq_out$burnout["late"])
    q1_C <- unname(seq_out$qC["first"])
    q2_C <- unname(seq_out$qC["late"])

    outcomes_by_strain <- c(
      both_burnout = Q1 * Q2,
      strain1_only_burnout = Q1 * (1 - Q2),
      strain2_only_burnout = (1 - Q1) * Q2,
      neither_burnout = (1 - Q1) * (1 - Q2)
    )

    seq_out$strain_burnout <- c(Q1 = Q1, Q2 = Q2)
    seq_out$qC_by_strain <- c(q1_C = q1_C, q2_C = q2_C)
    seq_out$G1_given_2 <- seq_out$G_first_given_late$G_at_qC
    seq_out$outcomes_by_strain <- outcomes_by_strain
    seq_out$boundary_event_found <- TRUE
    seq_out$boundary_event <- "yes"
    seq_out$entry_case <- "strain1_first_strain2_later"
    seq_out$first_strain <- "strain1"
    seq_out$event_time <- ode_event$time
    seq_out$event_state <- event_state
    seq_out$ystar <- c(ystar1 = ystar1, ystar2 = ystar2)
    seq_out$inputs <- inputs
    seq_out$ode <- ode_event
    return(seq_out)
  }

  if (first_strain == "strain2") {
    ## strain 2 enters first, strain 1 is late/resident. Compute by swapping labels.
    seq_out <- sequential_entry_wave(
      R0_first = R02,
      R0_late = R01,
      epsilon = epsilon,
      N = N,
      xa = x_event,
      y_late_a = y1_event,
      ystar_first = ystar2,
      ystar_late = ystar1,
      q_fun = q_fun,
      q_method = q_method,
      n_x = n_x,
      discrete_floor = FALSE
    )

    Q1 <- unname(seq_out$burnout["late"])
    Q2 <- unname(seq_out$burnout["first"])
    q1_C <- unname(seq_out$qC["late"])
    q2_C <- unname(seq_out$qC["first"])

    outcomes_by_strain <- c(
      both_burnout = Q1 * Q2,
      strain1_only_burnout = Q1 * (1 - Q2),
      strain2_only_burnout = (1 - Q1) * Q2,
      neither_burnout = (1 - Q1) * (1 - Q2)
    )

    seq_out$strain_burnout <- c(Q1 = Q1, Q2 = Q2)
    seq_out$qC_by_strain <- c(q1_C = q1_C, q2_C = q2_C)
    seq_out$G2_given_1 <- seq_out$G_first_given_late$G_at_qC
    seq_out$outcomes_by_strain <- outcomes_by_strain
    seq_out$boundary_event_found <- TRUE
    seq_out$boundary_event <- "yes"
    seq_out$entry_case <- "strain2_first_strain1_later"
    seq_out$first_strain <- "strain2"
    seq_out$event_time <- ode_event$time
    seq_out$event_state <- event_state
    seq_out$ystar <- c(ystar1 = ystar1, ystar2 = ystar2)
    seq_out$inputs <- inputs
    seq_out$ode <- ode_event
    return(seq_out)
  }

  stop("Unexpected event state.")
}
