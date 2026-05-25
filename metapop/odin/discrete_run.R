library(odin)
library(dde) ## odin insists on this
library(tidyr)
library(dplyr)
library(ggplot2)
library(future)
library(furrr)
library(patchwork)

s <- function(x) source(here::here("metapop/odin", x))
s("discrete_odin.R")
s("discrete_macpan2.R")
s("discrete_pureR.R")

##' @param K carrying capacity (scalar or vector of length n_patch)
##' @param r host growth rate per disease generation (ditto)
##' @param n_patch number of patches
##' @param n_strains number of strains [HARD-CODED]
##' @param nt time steps
##' @param alpha between-patch transmission  probability
##' @param I_init initial infection (chosen as Poisson random variable across platforms)
##' @param seed PRNG seed
##' @param chunk steps per chunk for early stopping (odin only; ignored if stop_cond is NULL)
##' @param stop_cond function(row) -> logical; called on the last row of each
##'   chunk; return TRUE to halt early. Use stop_either_extinct() or
##'   stop_both_extinct() from discrete_odin.R, or supply a custom function.
##'   NULL (default) runs all nt steps. Only supported for platform = "odin".
##' @param platform
##' @examples
##' ## stop as soon as either strain goes extinct (odin platform only):
##' ## plan(multisession)
##' ## runs <- discrete_run(nsim = 20, stop_cond = stop_either_extinct)
##' ## sumfun_discrete(runs)
##'
##' ## custom condition — stop when strain 2 alone is gone:
##' ## runs <- discrete_run(nsim = 20, chunk = 100,
##' ##   stop_cond = function(row) sum(row[, grep(",2\\]$", colnames(row))]) == 0)
## Parallelism is controlled by the caller via future::plan() before invoking
## discrete_run(). e.g. plan(multisession, workers = 8) for parallel runs,
## plan(sequential) (the default) for single-threaded execution.
discrete_run <- function(beta_vec = c(1.5, 2.5),
                          K = 1e4,
                          r = 0.125,
                          n_patch = 100,
                          nt = 1000,
                          alpha = 1e-3,
                          strain2_delay = 0,
                          I_init = 10,
                          seed = NULL,
                          nsim = 1,
                          chunk = NULL,
                          stop_cond = NULL,
                          platform = c("odin", "macpan2", "pureR")) {

  platform <- match.arg(platform)
  if ((!is.null(chunk) || !is.null(stop_cond)) && platform != "odin")
    stop("chunk and stop_cond are only supported for platform = 'odin'")

  args <- tibble::lst(beta_vec, r, K, n_patch, nt, I_init, alpha, strain2_delay)

  makefun  <- get(sprintf("make_simulator_%s", platform))
  runfun   <- get(sprintf("run_simulator_%s", platform))
  convfun  <- get(sprintf("conv_%s", platform))
  run_args <- Filter(Negate(is.null), list(chunk = chunk, stop_cond = stop_cond))

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
    select(-state)
}

## Analogue of sumfun in poisson/simulation_funs.R, adapted for the long-format
## list output of discrete_run() and extended to two strains.
##
## For each strain j:
##   - extinction time: first step where total Ij across patches == 0
##   - QE means: average of last `nsteps` steps over (step, run) cells where
##     Ij is still present, matching the cell-level masking of the original
##
## @param runs list of long-format tibbles from discrete_run() (or a single tibble)
## @param nsteps number of trailing steps used for quasi-equilibrium averages
## @return named numeric vector of summary statistics
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

## scratch work
if (FALSE) {
    
    ## no parallelization
    system.time(
        x <- discrete_run(nsim = 10)  ## 6.5 seconds
    )
    plan(multisession(workers = 4))
    system.time(
        x <- discrete_run(nsim = 10)  ## 1.23 seconds
    )
    sumfun_discrete(x)
  
}

