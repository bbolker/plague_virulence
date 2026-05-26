##' Run one or more two-strain n-patch simulations across platforms
##' @param beta_vec length-2 vector of per-capita transmission rates (one per strain)
##' @param K carrying capacity (scalar or vector of length n_patch)
##' @param r host growth rate per disease generation (ditto)
##' @param n_patch number of patches
##' @param nt time steps
##' @param alpha between-patch transmission probability
##' @param gamma length-2 vector of per-strain recovery rates (euler models only)
##' @param strain2_delay steps before strain 2 is seeded (0 = immediate)
##' @param I_init mean initial infected per patch (Poisson draw; length 1 or 2)
##' @param seed PRNG seed
##' @param nsim number of simulations to run
##' @param chunk steps per chunk for early stopping (odin only; ignored if stop_cond is NULL)
##' @param stop_cond function(row) -> logical; called on the last row of each
##'   chunk; return TRUE to halt early. Use stop_either_extinct() or
##'   stop_both_extinct(), or supply a custom function.
##'   NULL runs all nt steps. Only supported for platform = "odin".
##' @param def_file odin DSL filename in inst/odin/ (odin platform only)
##' @param dt time-step size in disease-generation units (euler models only)
##' @examples
##' ## stop as soon as either strain goes extinct (odin platform only):
##' ## plan(multisession)
##' ## runs <- discrete_run(nsim = 20, stop_cond = stop_either_extinct)
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
                         r = 0.125,
                         n_patch = 100,
                         nt = 1000,
                         alpha = 1e-3,
                         gamma = c(1, 1),
                          strain2_delay = 0,
                          I_init = 10,
                          seed = NULL,
                          nsim = 1,
                          chunk = NULL,
                          stop_cond = stop_both_extinct,
                          def_file = "discrete_odin_def.R",
                          dt = 1,
                          platform = c("odin", "macpan2", "pureR")) {

  platform <- match.arg(platform)
  if ((!is.null(chunk) || !is.null(stop_cond)) && platform != "odin")
    stop("chunk and stop_cond are only supported for platform = 'odin'")

  args <- tibble::lst(beta_vec, r, K, n_patch, nt, I_init, alpha, strain2_delay,
                      gamma, dt, def_file)

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
  FUN <- function(i) convfun(do.call(runfun, c(list(do.call(makefun, args)), run_args)))

  if (nsim == 1) {
    if (!is.null(seed)) set.seed(seed)
    return(FUN(1L))
  }

  furrr::future_map(seq.int(nsim), FUN,
                    .options = furrr::furrr_options(seed = seed %||% TRUE))
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
##' @param which strains to summarise (currently ignored; both strains always returned)
##' @return named numeric vector of summary statistics
##' @export
sumfun_discrete <- function(runs, nsteps = 100, which = c(1, 2)) {
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

  ## row index of first step where total infected == 0 (NA if never extinct)
  first_ext <- function(pop) {
    idx <- which(pop == 0)
    if (length(idx) == 0L) NA_integer_ else idx[1L]
  }
  ext1 <- sapply(agg, \(a) first_ext(a$I1_pop))
  ext2 <- sapply(agg, \(a) first_ext(a$I2_pop))

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
