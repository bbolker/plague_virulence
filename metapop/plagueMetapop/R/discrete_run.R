##' Construct the n_patch x n_strain matrix of initial infected counts
##' @description
##'   Expands \code{I_init} to a full \code{n_patch x n_strain} matrix.
##'   \code{I_init} may be a scalar (same value for all patches and strains),
##'   a length-\code{n_strain} vector (one value per strain, replicated across
##'   patches), or a full \code{n_patch x n_strain} matrix.
##' @param I_init scalar, length-\code{n_strain} vector, or \code{n_patch x n_strain} matrix
##' @param n_patch number of patches
##' @param n_strain number of strains (default 2)
##' @param method \code{"rpois"} (default): independent Poisson draws with
##'   \code{I_init} as per-cell means; \code{"fixed"}: use rounded \code{I_init}
##'   values directly (deterministic)
##' @return integer \code{n_patch x n_strain} matrix
##' @export
make_I_ini_mat <- function(I_init, n_patch, n_strain = 2L,
                            method = c("rpois", "fixed")) {
  method <- match.arg(method)
  if (is.matrix(I_init)) {
    stopifnot(nrow(I_init) == n_patch, ncol(I_init) == n_strain)
    lambda <- I_init
  } else {
    lambda <- matrix(rep(rep(I_init, length.out = n_strain), each = n_patch),
                     nrow = n_patch, ncol = n_strain)
  }
  if (method == "rpois")
    matrix(rpois(n_patch * n_strain, lambda = lambda), nrow = n_patch, ncol = n_strain)
  else {
    ## don't round lambda: up to user
    ## (warn if we can somehow figure out if a stoch model is being used and lambda is non-integer?)
    matrix(lambda, nrow = n_patch, ncol = n_strain)
  }
}

##' Run one or more two-strain n-patch simulations across platforms
##' @param beta_vec length-2 vector of per-capita transmission rates (one per strain)
##' @param K carrying capacity (scalar or vector of length n_patch)
##' @param logistic_growth 1 (default) = standard logistic demography; 0 = linear restoring force
##' @param reedfrost 1 = Reed-Frost mode (100\% removal per step); requires gamma == 1 and dt == 1
##' @param r host growth rate per disease generation (ditto)
##' @param n_patch number of patches
##' @param nt time steps
##' @param alpha between-patch transmission probability
##' @param gamma length-2 vector of per-strain recovery rates (euler models only)
##' @param strain2_delay steps before strain 2 is seeded (0 = immediate)
##' @param I_init scalar, length-2 vector, or n_patch x 2 matrix of initial infected
##'   counts (interpreted as means for \code{I_ini_method = "rpois"}, or used
##'   directly for \code{"fixed"}).  See \code{\link{make_I_ini_mat}}.
##' @param I_ini_method \code{"rpois"} (default): each simulation draws a fresh
##'   \code{n_patch x 2} matrix of Poisson counts; \code{"fixed"}: all simulations
##'   share the same rounded \code{I_init} values.
##' @param seed PRNG seed
##' @param nsim number of simulations to run
##' @param chunk steps per chunk for early stopping (odin only; ignored if stop_cond is NULL)
##' @param stop_cond function(row) -> logical; called on the last row of each
##'   chunk; return TRUE to halt early. Use stop_either_extinct() or
##'   stop_both_extinct(), or supply a custom function.
##'   NULL runs all nt steps. Only supported for platform = "odin".
##' @param def_file odin DSL filename in inst/odin/ (odin platform only).
##'   \code{"discrete_odin_def.R"}: stochastic discrete-time model (default).
##'   \code{"euler_odin_def.R"}: stochastic Euler model with explicit \code{dt},
##'     \code{gamma}, \code{logistic_growth}, \code{reedfrost}, and step-based
##'     strain-2 seeding via \code{strain2_delay}.
##'   \code{"euler_det_odin_def.R"}: deterministic Euler version (same step structure;
##'     flows replace random draws; still uses explicit \code{dt} and \code{strain2_delay}).
##'   \code{"ode_odin_def.R"}: deterministic ODE model using \code{deriv()} and odin's
##'     ODE solver; no \code{dt}, \code{strain2_delay}, or \code{reedfrost};
##'     \code{stop_cond} is not supported.
##' @param dt time-step size in disease-generation units (euler models only; ignored for ode_odin_def.R)
##' @examples
##' ## stop as soon as either strain goes extinct (odin platform only):
##' ## plan(multisession)
##' ## runs <- discrete_run(nsim = 20, stop_cond = stop_either_extinct())
##' ## sumfun_discrete(runs)
##'
##' ## custom condition — stop when strain 2 alone is gone:
##' ## runs <- discrete_run(nsim = 20, chunk = 100,
##' ##   stop_cond = function(row) sum(row[, grep(",2\\]$", colnames(row))]) == 0)
##' @param platform simulation backend: "odin" (default), "macpan2", or "pureR"
##' @return a single long-format tibble (nsim == 1) or a list of nsim tibbles
##' @details Parallelism is controlled by the caller via future::plan() before
##'   invoking discrete_run(). Use plan(multisession, workers = N) for parallel
##'   runs; plan(sequential) (the default) for single-threaded execution.
##' @export
discrete_run <- function(beta_vec = c(1.5, 2.5),
                         K = 1e4,
                         logistic_growth = 1,
                         reedfrost = 0,
                         r = if (logistic_growth) 0.125 else 0.04,
                         n_patch = 100,
                         nt = 1000,
                         alpha = 1e-3,
                         gamma = c(1, 1),
                          strain2_delay = 0,
                          I_init = 10,
                          I_ini_method = c("rpois", "fixed"),
                          seed = NULL,
                          nsim = 1,
                          chunk = NULL,
                          stop_cond = stop_both_extinct,
                          def_file = "discrete_odin_def.R",
                          dt = 1,
                          platform = c("odin", "macpan2", "pureR")) {

  .call        <- match.call()
  platform     <- match.arg(platform)
  I_ini_method <- match.arg(I_ini_method)
  .params      <- tibble::lst(beta_vec, K, logistic_growth, reedfrost, r, n_patch,
                               nt, alpha, gamma, strain2_delay, I_init, I_ini_method,
                               seed, nsim, chunk, stop_cond, def_file, dt, platform)
  attach_meta  <- function(x) {
    class(x) <- c("metapop_run", class(x))
    attr(x, "call")   <- .call
    attr(x, "params") <- .params
    x
  }
  if ((!is.null(chunk) || !is.null(stop_cond)) && platform != "odin")
    stop("chunk and stop_cond are only supported for platform = 'odin'")
  if (!is.null(stop_cond) && platform == "odin" && def_file == "ode_odin_def.R")
    stop("stop_cond is not supported for the deterministic ODE model (ode_odin_def.R)")
  if (reedfrost == 1 && !(all(gamma == 1) && dt == 1))
    stop("reedfrost = 1 requires gamma == 1 and dt == 1")

  args <- tibble::lst(beta_vec, r, K, n_patch, nt, alpha, strain2_delay,
                      gamma, dt, def_file)
  if (platform == "odin") {
    if (def_file %in% c("euler_odin_def.R", "ode_odin_def.R"))
      args$logistic_growth <- logistic_growth
    if (def_file == "euler_odin_def.R")
      args$reedfrost <- reedfrost
  }

  makefun  <- get(sprintf("make_simulator_%s", platform))
  runfun   <- get(sprintf("run_simulator_%s", platform))
  convfun  <- get(sprintf("conv_%s", platform))
  run_args <- Filter(Negate(is.null), list(chunk = chunk, stop_cond = stop_cond))

  ## pre-compile the odin generator once on the main process; workers receive it
  ## via the closure and skip recompilation (odin caches the .so by file hash,
  ## so even if serialization forces a reload, compilation is skipped)
  if (platform == "odin") args$gen_local <- compile_odin(def_file)

  ## Each future builds its own simulator. A per-worker cache (<<- into the
  ## closure) would be more efficient but furrr re-serializes the closure for
  ## every task, so the cache never persists across futures. If compilation
  ## time dominates, consider pre-building via makefun and passing the object
  ## as a future global (if odin objects serialize safely across workers).
  ## I_ini_mat is constructed inside FUN so each parallel simulation gets
  ## independent Poisson draws (for I_ini_method = "rpois"); for "fixed" all
  ## simulations share the same rounded values, but the copy is still made
  ## per-call to keep the interface uniform.
  FUN <- function(i) {
    local_args <- c(args,
                    list(I_ini_mat = make_I_ini_mat(I_init, n_patch, 2L, I_ini_method)))
    fun <- do.call(makefun, local_args)
    res <- do.call(runfun, c(list(fun), run_args)) |> convfun()
    if (dt != 1 && platform != "odin") dplyr::mutate(res, dplyr::across(step, ~ . * dt)) else res
  }

  if (nsim == 1) {
    if (!is.null(seed)) set.seed(seed)
    return(attach_meta(FUN(1L)))
  }

  attach_meta(furrr::future_map(seq.int(nsim), FUN,
                                .options = furrr::furrr_options(seed = seed %||% TRUE)))
}



##' Endemic equilibria (S* and I*) for the single-strain single-patch ODE model
##'
##' At the non-trivial equilibrium, dI/dt = 0 gives S* = gamma*K/beta.
##' Setting dS/dt = 0 then determines I*:
##' logistic: r*S*(1-S/K) = beta*I*S/K => I* = r*(K - S*)/beta;
##' linear:   r*(K-S)      = beta*I*S/K => I* = r*(K - S*)*K/(beta*S*).
##' Returns NA for both components when R0 = beta/gamma <= 1.
##' @param beta transmission rate (beta_vec[1])
##' @param gamma recovery rate (gamma[1]; default 1)
##' @param K carrying capacity
##' @param r host intrinsic growth rate
##' @param logistic_growth 1 = standard logistic (default); 0 = linear restoring force
##' @return named numeric vector with elements \code{eq_S} (S*) and \code{eq_I} (I1*),
##'   or \code{c(eq_S = NA, eq_I = NA)} when R0 <= 1
##' @export
ode_eq <- function(beta, gamma = 1, K = 1e4, r = 0.125, logistic_growth = 1) {
  if (beta <= gamma) return(c(eq_S = NA_real_, eq_I = NA_real_))
  S_star <- gamma * K / beta
  I_star <- if (logistic_growth == 1) {
    r * (K - S_star) / beta
  } else {
    r * (K - S_star) * K / (beta * S_star)
  }
  c(eq_S = S_star, eq_I = I_star)
}

##' Epidemic transient summary statistics for a single-strain single-patch ODE run
##'
##' @param run long-format tibble from \code{discrete_run()} with \code{n_patch = 1}
##'   and \code{I_init = c(I0, 0)} (single strain)
##' @param beta transmission rate (\code{beta_vec[1]}); taken from \code{attr(run, "params")} if NULL
##' @param gamma recovery rate (\code{gamma[1]}); taken from \code{attr(run, "params")} if NULL
##' @param K carrying capacity; taken from \code{attr(run, "params")} if NULL
##' @param r host intrinsic growth rate; taken from \code{attr(run, "params")} if NULL
##' @param logistic_growth 1 = logistic; 0 = linear restoring force; taken from \code{attr(run, "params")} if NULL
##' @return named numeric vector with elements:
##'   \code{eq} (endemic equilibrium of I1),
##'   \code{t_enter.boundary} (first downward crossing of eq by I1: last step with I1 > eq before dip),
##'   \code{t_Imin}, \code{Imin} (time and value of first local minimum of I1,
##'     identified by the first index where diff(diff(I1)) > 0 after \code{t_enter.boundary}),
##'   \code{t_Smin}, \code{Smin} (time and value of first local minimum of S),
##'   \code{t_eqS} (first time after \code{t_Smin} that S crosses upward through \code{eq_S}),
##'   \code{t_leave.boundary} (second upward crossing of eq by I1: last step with I1 < eq before recovery).
##'   \code{trough_area} (integral of \code{eq_S - S} from \code{t_enter.boundary}
##'     until the first time S rises back above \code{eq_S}).
##'   Any statistic that cannot be found returns \code{NA}.
##' @export
traj_stats_ode <- function(run, beta = NULL, gamma = NULL, K = NULL, r = NULL,
                           logistic_growth = NULL) {
  p <- attr(run, "params")
  if (is.null(beta))            beta            <- p$beta_vec[1L]
  if (is.null(gamma))           gamma           <- p$gamma[1L]
  if (is.null(K))               K               <- p$K[1L]
  if (is.null(r))               r               <- p$r[1L]
  if (is.null(logistic_growth)) logistic_growth <- p$logistic_growth
  agg <- dplyr::arrange(sum_run1(run, "wide"), step)
  I1  <- agg$I1_pop
  S   <- agg$S_pop
  t   <- agg$step
  n   <- length(t)

  eq   <- ode_eq(beta, gamma, K, r, logistic_growth)
  eq_S <- eq[["eq_S"]]
  eq_I <- eq[["eq_I"]]

  ## T1: first downward crossing (last step with I1 > eq_I before I1 drops below eq_I)
  down_idx <- which(I1[-n] > eq_I & I1[-1] < eq_I)
  T1 <- if (length(down_idx) > 0L) t[down_idx[1L]] else NA_real_

  ## T2, I2: first local minimum of I1 after T1
  cand2 <- which(diff(sign(diff(I1)))==2)
  if (!is.na(T1)) cand2 <- cand2[t[cand2] > T1]
  if (length(cand2) > 0L) {
    T2     <- t[cand2[1L]]
    I2_val <- I1[cand2[1L]]
  } else {
    T2 <- I2_val <- NA_real_
  }

  ## T3, I3: first local minimum of S
  cand3 <- which(diff(sign(diff(S)))==2)
  if (length(cand3) > 0L) {
    T3     <- t[cand3[1L]]
    I3_val <- S[cand3[1L]]
  } else {
    T3 <- I3_val <- NA_real_
  }

  ## T4: second upward crossing (last step with I1 < eq_I before I1 rises above eq_I)
  up_idx <- which(I1[-n] < eq_I & I1[-1] > eq_I)
  T4 <- if (length(up_idx) >= 2L) t[up_idx[2L]] else NA_real_

  ## T_S_recover: first upward crossing of eq_S by S (S recovers above equilibrium)
  S_up_idx    <- which(S[-n] < eq_S & S[-1] > eq_S)
  T_S_recover <- if (length(S_up_idx) > 0L) t[S_up_idx[1L]] else NA_real_

  ## t_eqS: first upward S crossing of eq_S after t_Smin
  S_up_after_Smin <- if (!is.na(T3)) S_up_idx[t[S_up_idx] > T3] else S_up_idx
  t_eqS <- if (length(S_up_after_Smin) > 0L) t[S_up_after_Smin[1L]] else NA_real_

  ## trough_area: integral of (eq_S - S) * dt from t_enter.boundary until S first exceeds eq_S
  if (!is.na(T1) && !is.na(T_S_recover)) {
    dt_step     <- if (n > 1L) t[2L] - t[1L] else 1
    trough_idx  <- which(t >= T1 & t <= T_S_recover)
    trough_area <- sum(eq_S - S[trough_idx]) * dt_step
  } else {
    trough_area <- NA_real_
  }

  c(eq_S = eq_S, eq_I = eq_I, t_enter.boundary = T1, t_Imin = T2, Imin = I2_val,
    t_Smin = T3, Smin = I3_val, t_eqS = t_eqS, t_leave.boundary = T4,
    trough_area = trough_area)
}



##' Aggregate a single run to per-step population totals
##' @param x long-format tibble from discrete_run()
##' @param return_type "long" or "wide"
##' @export
sum_run1 <- function(x, return_type = c("long", "wide")) {
  return_type <- match.arg(return_type)
  ret <- x |>
    dplyr::summarise(
      S_pop    = sum(value[state == "S"]),
      I1_pop   = sum(value[state == "I1"]),
      I2_pop   = sum(value[state == "I2"]),
      I1_patch = sum(value[state == "I1"] > 0),
      I2_patch = sum(value[state == "I2"] > 0),
      .by = step
    )
  if (return_type == "wide") return(ret)
  ret |>
    tidyr::pivot_longer(cols = -step, names_to = "state") |>
    dplyr::mutate(var  = sub(".*_", "", state),
                  type = sub("_.*", "", state)) |>
    dplyr::select(-state)
}

##' Summarise extinction and quasi-equilibrium statistics across runs
##' @details Analogue of sumfun in poisson/simulation_funs.R, adapted for the
##'   long-format list output of discrete_run() and extended to two strains.
##'   For each strain: extinction time (first step with total Ij == 0) and
##'   quasi-equilibrium means over the last \code{nsteps} steps, conditioned
##'   on the strain still being present (cell-level masking).
##' @param runs list of long-format tibbles from discrete_run(), or a single tibble
##' @param nsteps number of trailing steps used for quasi-equilibrium averages
##' @param strain2_delay steps before strain 2 is seeded; zeros in I2 before
##'   this step are not counted as extinction
##' @param which strains to summarise (currently ignored; both strains always returned)
##' @return named numeric vector of summary statistics
##' @export
sumfun_discrete <- function(runs, nsteps = 100, strain2_delay = 0, which = c(1, 2)) {
  if (!is.list(runs) || inherits(runs, "data.frame")) runs <- list(runs)
  nsim <- length(runs)

  ## aggregate each run to wide per-step format and sort by step
  agg <- lapply(runs, \(x) dplyr::arrange(sum_run1(x, "wide"), step))

  ## With early stopping, runs may have different lengths. Use the longest run
  ## to anchor the QE window. Shorter (extinct) runs return NA for out-of-bounds
  ## row accesses, but those runs are already masked by ext_vec, so the mean is
  ## unaffected.
  nt     <- max(vapply(agg, nrow, integer(1L)))
  window <- seq(max(1L, nt - nsteps + 1L), nt)

  ## row index of first step where total infected == 0 (NA if never extinct).
  ## `after`: ignore zeros at or before this step (used to skip pre-seeding
  ## zeros for strain 2).
  first_ext <- function(pop, steps = seq_along(pop), after = 0L) {
    idx <- which(pop == 0 & steps > after)
    if (length(idx) == 0L) NA_integer_ else idx[1L]
  }
  ext1 <- sapply(agg, \(a) first_ext(a$I1_pop))
  ext2 <- sapply(agg, \(a) first_ext(a$I2_pop, steps = a$step, after = strain2_delay))

  ## cell-level QE mean: mask (window-row, run) cells at or after extinction,
  ## then average over all remaining cells (matching original sumfun logic)
  qe_mean <- function(var, ext_vec) {
    m <- sapply(seq_len(nsim), \(i) {
      v <- agg[[i]][[var]][window]
      if (!is.na(ext_vec[i])) v[window >= ext_vec[i]] <- NA_real_
      v
    })
    mean(m, na.rm = TRUE)
  }

  c(
    mean_ext_time.I1 = mean(ext1, na.rm = TRUE),
    mean_ext_time.I2 = mean(ext2, na.rm = TRUE),
    ext_prob.I1      = mean(!is.na(ext1)),
    ext_prob.I2      = mean(!is.na(ext2)),
    ## QE means conditioned on I1 surviving
    qe_pop_S.I1  = qe_mean("S_pop",    ext1),
    qe_infpop.I1    = qe_mean("I1_pop",   ext1),
    qe_infpatch.I1  = qe_mean("I1_patch", ext1),
    ## QE means conditioned on I2 surviving
    qe_pop_S.I2  = qe_mean("S_pop",    ext2),
    qe_infpop.I2    = qe_mean("I2_pop",   ext2),
    qe_infpatch.I2  = qe_mean("I2_patch", ext2)
  )
}
