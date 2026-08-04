## Full stochastic validation utilities for the logistic burnout approximation.
## These functions deliberately leave the analytical approximation unchanged.

source(file.path("fadeout", "logistic_burnout", "logistic_burnout_functions.R"))

.require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required")
  }
}

parse_cli_options <- function(defaults = list()) {
  args <- commandArgs(trailingOnly = TRUE)
  out <- defaults
  if (!length(args)) return(out)
  for (arg in args) {
    if (!grepl("^--", arg)) next
    z <- sub("^--", "", arg)
    key <- sub("=.*$", "", z)
    val <- if (grepl("=", z)) sub("^[^=]*=", "", z) else "TRUE"
    out[[key]] <- val
  }
  out
}

as_num <- function(x) as.numeric(x)
as_int <- function(x) as.integer(as.numeric(x))
as_logical <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

format_num <- function(x, digits = 4) {
  ifelse(is.finite(x), format(signif(x, digits), scientific = FALSE), NA)
}

wilson_ci <- function(x, n, conf.level = 0.95) {
  if (is.na(x) || is.na(n) || n <= 0) {
    return(c(estimate = NA_real_, se = NA_real_, low = NA_real_, high = NA_real_))
  }
  p <- x / n
  z <- stats::qnorm(1 - (1 - conf.level) / 2)
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n) / denom
  c(
    estimate = p,
    se = sqrt(p * (1 - p) / n),
    low = max(0, centre - half),
    high = min(1, centre + half)
  )
}

ci_half_width <- function(x, n) {
  ci <- wilson_ci(x, n)
  (ci[["high"]] - ci[["low"]]) / 2
}

tau_delta <- function(R0, I0 = 1, delta = 1e-6) {
  if (!is.finite(R0) || R0 <= 1) return(NA_real_)
  if (!is.finite(I0) || I0 < 1) stop("I0 must be >= 1")
  if (!is.finite(delta) || delta <= 0 || delta >= 1) {
    stop("delta must be in (0, 1)")
  }
  a <- exp(-log1p(-delta) / I0)
  numerator <- a - 1 / R0
  denominator <- a - 1
  if (numerator <= 0 || denominator <= 0) return(NA_real_)
  log(numerator / denominator) / (R0 - 1)
}

analytical_validation_quantities <- function(R0, r, K, I0 = 1) {
  ans <- logistic_burnout_probability(
    R0 = R0, r = r, K = K, I0 = I0,
    lineage_count_method = "round",
    initial_tmax = 50, maximum_tmax = 1600, dt = 0.02,
    rtol = 1e-9, atol = 1e-11,
    rel.tol = 1e-8, subdivisions = 500
  )
  Q <- ans$P_burnout
  P_cond <- if (is.finite(Q)) 1 - Q else NA_real_
  p_est <- 1 - (1 / R0)^I0
  P_uncond <- p_est * P_cond
  m_used <- ans$m_used
  if (!is.finite(m_used)) {
    m_used <- max(1, round(ans$m_raw))
  }
  data.frame(
    x_in = ans$x_in,
    y_star = ans$y_star,
    m_raw = ans$m_raw,
    m_used = m_used,
    q1 = ans$q1,
    Q_approx = Q,
    P_cond_approx = P_cond,
    p_est_approx = p_est,
    P_uncond_approx = P_uncond,
    analytical_status = ans$status,
    boundary_entry_found = ans$boundary_entry_found,
    integration_converged = ans$integration_converged,
    stringsAsFactors = FALSE
  )
}

make_parameter_grid <- function(mode) {
  if (mode == "smoke") {
    return(expand.grid(
      R0 = c(1.1, 2.5, 5),
      r = 0.1,
      K = c(1000, 10000),
      I0 = 1,
      stringsAsFactors = FALSE
    ))
  }
  if (mode == "main") {
    R0_log <- 1 + exp(seq(log(0.02), log(4), length.out = 24))
    R0_values <- sort(unique(round(c(
      R0_log, 1.05, 1.1, 1.5, 2, 2.5, 3, 4, 5
    ), 5)))
    return(expand.grid(
      R0 = R0_values,
      r = 0.1,
      K = c(1000, 3000, 10000, 30000),
      I0 = 1,
      stringsAsFactors = FALSE
    ))
  }
  if (mode == "r_sensitivity") {
    return(expand.grid(
      R0 = c(1.05, 1.1, 1.5, 2, 2.5, 3, 4, 5),
      r = c(0.05, 0.1, 0.125, 0.2),
      K = c(3000, 10000),
      I0 = 1,
      stringsAsFactors = FALSE
    ))
  }
  if (mode == "exact_subset") {
    return(data.frame(
      R0 = c(1.05, 1.5, 2.5, 5, 2.5),
      r = c(0.1, 0.1, 0.1, 0.1, 0.1),
      K = c(10000, 10000, 10000, 1000, 3000),
      I0 = 1,
      stringsAsFactors = FALSE
    ))
  }
  stop("Unknown mode: ", mode)
}

cell_seed <- function(global_seed, R0, r, K, I0, extra = 0L) {
  raw <- as.integer(round(c(R0 * 100000, r * 100000, K, I0, extra)))
  as.integer((as.numeric(global_seed) + sum(raw * c(3, 5, 7, 11, 13))) %% 
               .Machine$integer.max)
}

simulate_tau_leap_full <- function(R0, r, K, I0 = 1, n_sim = 1000,
                                   dt = 0.02, delta = 1e-6,
                                   seed = 1L, tmax = 200) {
  set.seed(seed)
  td <- tau_delta(R0, I0, delta)
  if (!is.finite(td)) stop("tau_delta is not finite")
  analytical <- analytical_validation_quantities(R0, r, K, I0)
  I_BL <- as.integer(analytical$m_used[1])
  if (!is.finite(I_BL) || I_BL < 1) stop("Invalid analytical boundary count")

  S <- rep(as.integer(round(K - I0)), n_sim)
  I <- rep(as.integer(round(I0)), n_sim)
  status <- rep("active", n_sim)
  phase <- rep("early", n_sim)
  above_boundary <- rep(FALSE, n_sim)
  entered_boundary <- rep(FALSE, n_sim)
  below_after_entry <- rep(FALSE, n_sim)
  resolved_time <- rep(NA_real_, n_sim)
  p_rec <- -expm1(-dt)
  steps <- ceiling(tmax / dt)

  for (step in seq_len(steps)) {
    active <- which(status == "active")
    if (!length(active)) break
    time <- step * dt
    prev_I <- I[active]

    extinct_pre <- active[I[active] <= 0L]
    if (length(extinct_pre)) {
      early_ext <- extinct_pre[phase[extinct_pre] == "early" & time <= td]
      late_ext <- setdiff(extinct_pre, early_ext)
      status[early_ext] <- "fizzle"
      status[late_ext] <- "late_extinction_before_boundary"
      resolved_time[c(early_ext, late_ext)] <- time
      active <- which(status == "active")
      if (!length(active)) break
      prev_I <- I[active]
    }

    hazard <- R0 * I[active] / K
    p_inf <- -expm1(-hazard * dt)
    incidence <- stats::rbinom(length(active), S[active], pmin(1, pmax(0, p_inf)))
    removals <- stats::rbinom(length(active), I[active], p_rec)
    delta_log <- r * S[active] * (1 - S[active] / K) * dt
    pop_change <- integer(length(active))
    pos <- delta_log >= 0
    pop_change[pos] <- stats::rpois(sum(pos), pmax(0, delta_log[pos]))
    neg <- !pos
    if (any(neg)) {
      prob <- pmin(1, pmax(0, -delta_log[neg] / pmax(1, S[active][neg])))
      pop_change[neg] <- -stats::rbinom(sum(neg), S[active][neg], prob)
    }
    S[active] <- pmax(0L, S[active] - incidence + pop_change)
    I[active] <- pmax(0L, I[active] + incidence - removals)

    if (time >= td) phase[active[phase[active] == "early"]] <- "established"
    established <- active[phase[active] != "early"]
    above_boundary[established] <- above_boundary[established] |
      I[established] > I_BL

    crossed_down <- established[
      above_boundary[established] &
        !entered_boundary[established] &
        prev_I[match(established, active)] > I_BL &
        I[established] <= I_BL
    ]
    if (length(crossed_down)) {
      entered_boundary[crossed_down] <- TRUE
      below_after_entry[crossed_down] <- TRUE
    }
    in_boundary <- established[entered_boundary[established]]
    below_after_entry[in_boundary] <- below_after_entry[in_boundary] |
      I[in_boundary] <= I_BL

    extinct <- active[I[active] <= 0L]
    if (length(extinct)) {
      fizzle <- extinct[phase[extinct] == "early" & time <= td]
      burnout <- extinct[entered_boundary[extinct]]
      late <- setdiff(extinct, c(fizzle, burnout))
      status[fizzle] <- "fizzle"
      status[burnout] <- "burnout"
      status[late] <- "late_extinction_before_boundary"
      resolved_time[c(fizzle, burnout, late)] <- time
    }
    candidates <- in_boundary[
      status[in_boundary] == "active" &
        below_after_entry[in_boundary] &
        I[in_boundary] > I_BL
    ]
    if (length(candidates)) {
      status[candidates] <- "persistence"
      resolved_time[candidates] <- time
    }
    if (any(S < 0) || any(I < 0)) stop("Negative state detected")
  }
  status[status == "active"] <- "unresolved"
  make_full_summary(status, n_sim, R0, r, K, I0, delta, td, seed,
                    engine = "tau_leap", dt = dt, analytical = analytical,
                    tmax = tmax)
}

classify_extinction <- function(phase, time, tau_delta_value, entered_boundary) {
  if (phase == "early" && time <= tau_delta_value) return("fizzle")
  if (entered_boundary) return("burnout")
  "late_extinction_before_boundary"
}

simulate_adaptive_tau_full_one <- function(
    R0, r, K, I0 = 1, delta = 1e-6, seed = NULL, tmax = 200,
    rel_tol = 0.03, max_tau = 0.25, min_tau = 1e-5,
    exact_I_threshold = 100L, exact_rate_threshold = 250) {
  if (!is.null(seed)) set.seed(seed)
  td <- tau_delta(R0, I0, delta)
  analytical <- analytical_validation_quantities(R0, r, K, I0)
  I_BL <- as.integer(analytical$m_used[1])
  S <- as.integer(round(K - I0))
  I <- as.integer(round(I0))
  time <- 0
  phase <- "early"
  above <- FALSE
  entered <- FALSE
  below <- FALSE

  while (time < tmax) {
    if (I <= 0L) return(classify_extinction(phase, time, td, entered))
    if (time >= td && phase == "early") phase <- "established"

    infection_rate <- R0 * S * I / K
    removal_rate <- I
    growth_rate <- r * S * (1 - S / K)
    birth_rate <- max(0, growth_rate)
    sus_death_rate <- max(0, -growth_rate)
    total_rate <- infection_rate + removal_rate + birth_rate + sus_death_rate
    if (!is.finite(total_rate) || total_rate <= 0) break

    use_exact <- I <= max(exact_I_threshold, 2L * I_BL) ||
      total_rate <= exact_rate_threshold

    prev_I <- I
    if (use_exact) {
      time <- time + stats::rexp(1, total_rate)
      if (time > tmax) break
      u <- stats::runif(1) * total_rate
      if (u <= infection_rate) {
        if (S > 0) {
          S <- S - 1L
          I <- I + 1L
        }
      } else if (u <= infection_rate + removal_rate) {
        I <- max(0L, I - 1L)
      } else if (u <= infection_rate + removal_rate + birth_rate) {
        S <- S + 1L
      } else {
        S <- max(0L, S - 1L)
      }
    } else {
      tau_s <- rel_tol * max(S, 1) /
        max(infection_rate + birth_rate + sus_death_rate, .Machine$double.eps)
      tau_i <- rel_tol * max(I, 1) /
        max(infection_rate + removal_rate, .Machine$double.eps)
      tau <- min(max_tau, tau_s, tau_i, tmax - time)
      if (!is.finite(tau) || tau < min_tau) {
        tau <- min(min_tau, tmax - time)
      }
      inf_events <- stats::rpois(1, infection_rate * tau)
      rem_events <- stats::rpois(1, removal_rate * tau)
      birth_events <- stats::rpois(1, birth_rate * tau)
      sus_death_events <- stats::rpois(1, sus_death_rate * tau)
      inf_events <- min(inf_events, S)
      rem_events <- min(rem_events, I)
      sus_death_events <- min(sus_death_events, max(0L, S - inf_events))
      S <- S - inf_events + birth_events - sus_death_events
      I <- I + inf_events - rem_events
      S <- max(0L, as.integer(S))
      I <- max(0L, as.integer(I))
      time <- time + tau
    }

    if (phase != "early") {
      above <- above || I > I_BL
      if (above && !entered && prev_I > I_BL && I <= I_BL) {
        entered <- TRUE
        below <- TRUE
      }
      if (entered) below <- below || I <= I_BL
      if (entered && below && I > I_BL) return("persistence")
    }
    if (S < 0 || I < 0) stop("Negative state detected")
  }
  "unresolved"
}

simulate_adaptive_tau_full <- function(
    R0, r, K, I0 = 1, n_sim = 1000, delta = 1e-6,
    seed = 1L, tmax = 200, rel_tol = 0.03, max_tau = 0.25,
    boundary_tau = 0.01, exact_I_threshold = 100L,
    exact_rate_threshold = 250) {
  set.seed(seed)
  td <- tau_delta(R0, I0, delta)
  analytical <- analytical_validation_quantities(R0, r, K, I0)
  I_BL <- as.integer(analytical$m_used[1])
  if (!is.finite(I_BL) || I_BL < 1) stop("Invalid analytical boundary count")

  S <- rep(as.integer(round(K - I0)), n_sim)
  I <- rep(as.integer(round(I0)), n_sim)
  time <- rep(0, n_sim)
  status <- rep("active", n_sim)
  phase <- rep("early", n_sim)
  above_boundary <- rep(FALSE, n_sim)
  entered_boundary <- rep(FALSE, n_sim)
  below_after_entry <- rep(FALSE, n_sim)

  while (any(status == "active" & time < tmax)) {
    active <- which(status == "active" & time < tmax)
    phase[active[phase[active] == "early" & time[active] >= td]] <- "established"

    infection_rate <- R0 * S[active] * I[active] / K
    removal_rate <- I[active]
    growth_rate <- r * S[active] * (1 - S[active] / K)
    birth_rate <- pmax(0, growth_rate)
    sus_death_rate <- pmax(0, -growth_rate)
    total_rate <- infection_rate + removal_rate + birth_rate + sus_death_rate

    dead_or_stuck <- active[I[active] <= 0L | !is.finite(total_rate) |
                              total_rate <= 0]
    if (length(dead_or_stuck)) {
      for (idx in dead_or_stuck) {
        status[idx] <- classify_extinction(
          phase[idx], time[idx], td, entered_boundary[idx]
        )
      }
      active <- setdiff(active, dead_or_stuck)
      if (!length(active)) next
      infection_rate <- R0 * S[active] * I[active] / K
      removal_rate <- I[active]
      growth_rate <- r * S[active] * (1 - S[active] / K)
      birth_rate <- pmax(0, growth_rate)
      sus_death_rate <- pmax(0, -growth_rate)
      total_rate <- infection_rate + removal_rate + birth_rate + sus_death_rate
    }

    near_boundary <- phase[active] != "early" &
      (entered_boundary[active] | I[active] <= max(exact_I_threshold, 2L * I_BL))
    tau_s <- rel_tol * pmax(S[active], 1) /
      pmax(infection_rate + birth_rate + sus_death_rate, .Machine$double.eps)
    tau_i <- rel_tol * pmax(I[active], 1) /
      pmax(infection_rate + removal_rate, .Machine$double.eps)
    tau <- pmin(max_tau, tau_s, tau_i, tmax - time[active])
    tau[near_boundary] <- pmin(tau[near_boundary], boundary_tau)
    crosses_tau_delta <- phase[active] == "early" &
      time[active] < td & time[active] + tau > td
    tau[crosses_tau_delta] <- td - time[active[crosses_tau_delta]]
    tau <- pmax(tau, 1e-8)

    prev_I <- I[active]
    p_inf <- pmin(1, pmax(0, -expm1(-R0 * I[active] * tau / K)))
    p_rec <- pmin(1, pmax(0, -expm1(-tau)))
    incidence <- stats::rbinom(length(active), S[active], p_inf)
    removals <- stats::rbinom(length(active), I[active], p_rec)
    pop_change <- integer(length(active))
    pos <- growth_rate >= 0
    pop_change[pos] <- stats::rpois(sum(pos), birth_rate[pos] * tau[pos])
    neg <- !pos
    if (any(neg)) {
      prob <- pmin(1, pmax(
        0, sus_death_rate[neg] * tau[neg] / pmax(1, S[active][neg])
      ))
      pop_change[neg] <- -stats::rbinom(sum(neg), S[active][neg], prob)
    }
    S[active] <- pmax(0L, S[active] - incidence + pop_change)
    I[active] <- pmax(0L, I[active] + incidence - removals)
    time[active] <- time[active] + tau

    phase[active[phase[active] == "early" & time[active] >= td]] <- "established"
    established <- active[phase[active] != "early" & status[active] == "active"]
    above_boundary[established] <- above_boundary[established] |
      I[established] > I_BL

    prev_established_I <- prev_I[match(established, active)]
    crossed_down <- established[
      above_boundary[established] &
        !entered_boundary[established] &
        prev_established_I > I_BL &
        I[established] <= I_BL
    ]
    if (length(crossed_down)) {
      entered_boundary[crossed_down] <- TRUE
      below_after_entry[crossed_down] <- TRUE
    }
    in_boundary <- established[entered_boundary[established]]
    below_after_entry[in_boundary] <- below_after_entry[in_boundary] |
      I[in_boundary] <= I_BL

    extinct <- active[I[active] <= 0L & status[active] == "active"]
    if (length(extinct)) {
      for (idx in extinct) {
        status[idx] <- classify_extinction(
          phase[idx], time[idx], td, entered_boundary[idx]
        )
      }
    }
    candidates <- in_boundary[
      status[in_boundary] == "active" &
        below_after_entry[in_boundary] &
        I[in_boundary] > I_BL
    ]
    if (length(candidates)) status[candidates] <- "persistence"
    if (any(S < 0) || any(I < 0)) stop("Negative state detected")
  }
  status[status == "active"] <- "unresolved"
  make_full_summary(status, n_sim, R0, r, K, I0, delta,
                    td, seed,
                    engine = "adaptive_tau", dt = NA_real_,
                    analytical = analytical, tmax = tmax)
}

simulate_tau_leap_boundary <- function(R0, r, K, I0 = 1, n_sim = 1000,
                                       dt = 0.02, seed = 1L, tmax = 200) {
  set.seed(seed)
  analytical <- analytical_validation_quantities(R0, r, K, I0)
  if (!isTRUE(analytical$boundary_entry_found[1])) {
    return(make_boundary_summary(
      rep("unresolved", n_sim), n_sim, R0, r, K, I0, seed,
      engine = "tau_leap", dt = dt, analytical = analytical, tmax = tmax
    ))
  }
  I_BL <- as.integer(analytical$m_used[1])
  S <- rep(max(0L, as.integer(round(K * analytical$x_in[1]))), n_sim)
  I <- rep(I_BL, n_sim)
  status <- rep("active", n_sim)
  below_seen <- rep(TRUE, n_sim)
  p_rec <- -expm1(-dt)
  steps <- ceiling(tmax / dt)
  for (step in seq_len(steps)) {
    active <- which(status == "active")
    if (!length(active)) break
    hazard <- R0 * I[active] / K
    incidence <- stats::rbinom(length(active), S[active],
                               pmin(1, pmax(0, -expm1(-hazard * dt))))
    removals <- stats::rbinom(length(active), I[active], p_rec)
    delta_log <- r * S[active] * (1 - S[active] / K) * dt
    pop_change <- integer(length(active))
    pos <- delta_log >= 0
    pop_change[pos] <- stats::rpois(sum(pos), pmax(0, delta_log[pos]))
    neg <- !pos
    if (any(neg)) {
      prob <- pmin(1, pmax(0, -delta_log[neg] / pmax(1, S[active][neg])))
      pop_change[neg] <- -stats::rbinom(sum(neg), S[active][neg], prob)
    }
    S[active] <- pmax(0L, S[active] - incidence + pop_change)
    I[active] <- pmax(0L, I[active] + incidence - removals)
    below_seen[active] <- below_seen[active] | I[active] <= I_BL
    extinct <- active[I[active] <= 0L]
    exit <- active[I[active] > I_BL & below_seen[active]]
    status[extinct] <- "extinction"
    status[exit] <- "boundary_exit"
    if (any(S < 0) || any(I < 0)) stop("Negative state detected")
  }
  status[status == "active"] <- "unresolved"
  make_boundary_summary(status, n_sim, R0, r, K, I0, seed,
                        engine = "tau_leap", dt = dt,
                        analytical = analytical, tmax = tmax)
}

simulate_gillespie_full_one <- function(R0, r, K, I0 = 1, delta = 1e-6,
                                        seed = NULL, tmax = 200) {
  if (!is.null(seed)) set.seed(seed)
  td <- tau_delta(R0, I0, delta)
  analytical <- analytical_validation_quantities(R0, r, K, I0)
  I_BL <- as.integer(analytical$m_used[1])
  S <- as.integer(round(K - I0))
  I <- as.integer(round(I0))
  time <- 0
  phase <- "early"
  above <- FALSE
  entered <- FALSE
  below <- FALSE
  while (time < tmax) {
    if (I <= 0L) {
      if (phase == "early" && time <= td) return("fizzle")
      if (entered) return("burnout")
      return("late_extinction_before_boundary")
    }
    if (time >= td && phase == "early") phase <- "established"
    infection_rate <- R0 * S * I / K
    removal_rate <- I
    growth_rate <- r * S * (1 - S / K)
    birth_rate <- max(0, growth_rate)
    sus_death_rate <- max(0, -growth_rate)
    total_rate <- infection_rate + removal_rate + birth_rate + sus_death_rate
    if (!is.finite(total_rate) || total_rate <= 0) break
    time <- time + stats::rexp(1, total_rate)
    u <- stats::runif(1) * total_rate
    prev_I <- I
    if (u <= infection_rate) {
      if (S > 0) {
        S <- S - 1L
        I <- I + 1L
      }
    } else if (u <= infection_rate + removal_rate) {
      I <- max(0L, I - 1L)
    } else if (u <= infection_rate + removal_rate + birth_rate) {
      S <- S + 1L
    } else {
      S <- max(0L, S - 1L)
    }
    if (phase != "early") {
      above <- above || I > I_BL
      if (above && !entered && prev_I > I_BL && I <= I_BL) {
        entered <- TRUE
        below <- TRUE
      }
      if (entered) below <- below || I <= I_BL
      if (entered && below && I > I_BL) return("persistence")
    }
    if (S < 0 || I < 0) stop("Negative state detected")
  }
  "unresolved"
}

simulate_gillespie_full <- function(R0, r, K, I0 = 1, n_sim = 500,
                                    delta = 1e-6, seed = 1L, tmax = 200) {
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, n_sim)
  status <- vapply(seeds, function(s) {
    simulate_gillespie_full_one(R0, r, K, I0, delta, s, tmax)
  }, character(1))
  analytical <- analytical_validation_quantities(R0, r, K, I0)
  make_full_summary(status, n_sim, R0, r, K, I0, delta,
                    tau_delta(R0, I0, delta), seed,
                    engine = "gillespie", dt = NA_real_,
                    analytical = analytical, tmax = tmax)
}

make_full_summary <- function(status, n_total, R0, r, K, I0, delta, td, seed,
                              engine, dt, analytical, tmax) {
  counts <- table(factor(status, levels = c(
    "fizzle", "late_extinction_before_boundary", "burnout",
    "persistence", "unresolved"
  )))
  n_fizzle <- as.integer(counts[["fizzle"]])
  n_late <- as.integer(counts[["late_extinction_before_boundary"]])
  n_burnout <- as.integer(counts[["burnout"]])
  n_persistence <- as.integer(counts[["persistence"]])
  n_unresolved <- as.integer(counts[["unresolved"]])
  n_established <- n_total - n_fizzle
  n_boundary <- n_persistence + n_burnout
  ci_uncond <- wilson_ci(n_persistence, n_total)
  ci_fizzle <- wilson_ci(n_fizzle, n_total)
  ci_cond_est <- wilson_ci(n_persistence, n_established)
  ci_cond_bound <- wilson_ci(n_persistence, n_boundary)
  ci_q_bound <- wilson_ci(n_burnout, n_boundary)
  data.frame(
    R0 = R0, r = r, K = K, I0 = I0, delta = delta,
    tau_delta = td, engine = engine, dt = dt, seed = seed,
    tmax = tmax, n_total = n_total,
    n_fizzle = n_fizzle,
    n_late_extinction_before_boundary = n_late,
    n_burnout = n_burnout,
    n_persistence = n_persistence,
    n_unresolved = n_unresolved,
    p_fizzle_sim = ci_fizzle[["estimate"]],
    p_fizzle_sim_ci_low = ci_fizzle[["low"]],
    p_fizzle_sim_ci_high = ci_fizzle[["high"]],
    p_est_sim = 1 - ci_fizzle[["estimate"]],
    P_uncond_sim = ci_uncond[["estimate"]],
    P_uncond_sim_se = ci_uncond[["se"]],
    P_uncond_sim_ci_low = ci_uncond[["low"]],
    P_uncond_sim_ci_high = ci_uncond[["high"]],
    P_cond_sim_established = ci_cond_est[["estimate"]],
    P_cond_sim_established_ci_low = ci_cond_est[["low"]],
    P_cond_sim_established_ci_high = ci_cond_est[["high"]],
    P_cond_sim_boundary = ci_cond_bound[["estimate"]],
    P_cond_sim_boundary_ci_low = ci_cond_bound[["low"]],
    P_cond_sim_boundary_ci_high = ci_cond_bound[["high"]],
    Q_sim_boundary = ci_q_bound[["estimate"]],
    Q_sim_boundary_ci_low = ci_q_bound[["low"]],
    Q_sim_boundary_ci_high = ci_q_bound[["high"]],
    x_in = analytical$x_in,
    y_star = analytical$y_star,
    m_raw = analytical$m_raw,
    m_used = analytical$m_used,
    q1 = analytical$q1,
    Q_approx = analytical$Q_approx,
    P_cond_approx = analytical$P_cond_approx,
    p_est_approx = analytical$p_est_approx,
    P_uncond_approx = analytical$P_uncond_approx,
    error_unconditional = ci_uncond[["estimate"]] - analytical$P_uncond_approx,
    abs_error_unconditional = abs(ci_uncond[["estimate"]] - analytical$P_uncond_approx),
    error_conditional_established = ci_cond_est[["estimate"]] - analytical$P_cond_approx,
    error_conditional_boundary = ci_cond_bound[["estimate"]] - analytical$P_cond_approx,
    late_extinction_fraction = n_late / n_total,
    unresolved_fraction = n_unresolved / n_total,
    analytical_status = analytical$analytical_status,
    simulation_status = if (n_unresolved > 0) "has_unresolved" else "resolved",
    validity_notes = "",
    stringsAsFactors = FALSE
  )
}

make_boundary_summary <- function(status, n_total, R0, r, K, I0, seed,
                                  engine, dt, analytical, tmax) {
  counts <- table(factor(status, levels = c(
    "extinction", "boundary_exit", "unresolved"
  )))
  n_ext <- as.integer(counts[["extinction"]])
  n_exit <- as.integer(counts[["boundary_exit"]])
  n_unres <- as.integer(counts[["unresolved"]])
  denom <- n_ext + n_exit
  ci_q <- wilson_ci(n_ext, denom)
  data.frame(
    R0 = R0, r = r, K = K, I0 = I0,
    x_in = analytical$x_in,
    y_star = analytical$y_star,
    m_raw = analytical$m_raw,
    m_used = analytical$m_used,
    n_sim = n_total,
    n_extinction = n_ext,
    n_boundary_exit = n_exit,
    n_unresolved = n_unres,
    Q_boundary_start_sim = ci_q[["estimate"]],
    Q_boundary_start_se = ci_q[["se"]],
    Q_boundary_start_ci_low = ci_q[["low"]],
    Q_boundary_start_ci_high = ci_q[["high"]],
    Q_approx = analytical$Q_approx,
    error = ci_q[["estimate"]] - analytical$Q_approx,
    abs_error = abs(ci_q[["estimate"]] - analytical$Q_approx),
    engine = engine, dt = dt, seed = seed, tmax = tmax,
    analytical_status = analytical$analytical_status,
    stringsAsFactors = FALSE
  )
}

adaptive_full_cell <- function(R0, r, K, I0, delta, engine = "tau_leap",
                               dt = 0.02, seed = 1L,
                               initial_batch = 1000L,
                               min_sim = 5000L, max_sim = 50000L,
                               target_half_width = 0.01,
                               tmax = 200) {
  statuses <- character()
  total <- 0L
  repeat {
    batch <- min(initial_batch, max_sim - total)
    if (batch <= 0) break
    ans <- if (engine == "tau_leap") {
      simulate_tau_leap_full(R0, r, K, I0, batch, dt, delta,
                             seed + total + 1L, tmax)
    } else if (engine == "adaptive_tau") {
      simulate_adaptive_tau_full(R0, r, K, I0, batch, delta,
                                 seed + total + 1L, tmax)
    } else {
      simulate_gillespie_full(R0, r, K, I0, batch, delta,
                              seed + total + 1L, tmax)
    }
    batch_status <- rep("fizzle", ans$n_fizzle)
    batch_status <- c(batch_status,
      rep("late_extinction_before_boundary",
          ans$n_late_extinction_before_boundary),
      rep("burnout", ans$n_burnout),
      rep("persistence", ans$n_persistence),
      rep("unresolved", ans$n_unresolved))
    statuses <- c(statuses, batch_status)
    total <- length(statuses)
    n_persist <- sum(statuses == "persistence")
    if (total >= min_sim && ci_half_width(n_persist, total) <= target_half_width) {
      break
    }
    if (total >= max_sim) break
  }
  analytical <- analytical_validation_quantities(R0, r, K, I0)
  make_full_summary(
    statuses, length(statuses), R0, r, K, I0, delta,
    tau_delta(R0, I0, delta), seed, engine, dt, analytical, tmax
  )
}

check_full_table <- function(tab) {
  prob_cols <- names(tab)[grepl("(^p_|^P_|^Q_).*sim|approx|ci_", names(tab))]
  for (col in prob_cols) {
    x <- tab[[col]]
    bad <- is.finite(x) & (x < -1e-12 | x > 1 + 1e-12)
    if (any(bad)) stop("Probability outside [0,1] in ", col)
  }
  count_sum <- tab$n_fizzle + tab$n_late_extinction_before_boundary +
    tab$n_burnout + tab$n_persistence + tab$n_unresolved
  if (any(count_sum != tab$n_total)) stop("Outcome counts do not sum to n_total")
  if (any(abs(tab$P_cond_approx - (1 - tab$Q_approx)) > 1e-10, na.rm = TRUE)) {
    stop("Approximate conditional identity failed")
  }
  if (any(abs(tab$P_uncond_approx -
              tab$p_est_approx * tab$P_cond_approx) > 1e-10, na.rm = TRUE)) {
    stop("Approximate unconditional identity failed")
  }
  invisible(TRUE)
}

session_metadata <- function() {
  sha <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_
  )
  data.frame(
    git_sha = if (length(sha)) sha[1] else NA_character_,
    R_version = paste(R.version$major, R.version$minor, sep = "."),
    platform = R.version$platform,
    stringsAsFactors = FALSE
  )
}
