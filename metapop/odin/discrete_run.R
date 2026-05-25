library(odin)
library(dde) ## odin insists on this
library(tidyr)
library(dplyr)
library(ggplot2)
library(parallel)
library(patchwork)

source("discrete_odin.R")
source("discrete_macpan2.R")
source("discrete_pureR.R")

##' @param K carrying capacity (scalar or vector of length n_patch)
##' @param r host growth rate per disease generation (ditto)
##' @param n_patch number of patches
##' @param n_strains number of strains [HARD-CODED]
##' @param nt time steps
##' @param alpha between-patch transmission  probability
##' @param I_init initial infection (chosen as Poisson random variable across platforms)
##' @param seed PRNG seed
##' @param platform
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
                          cl = NULL,
                          ncores = 1,
                          platform = c("odin", "macpan2", "pureR")) {

  platform <- match.arg(platform)

  args <- tibble::lst(beta_vec, r, K, n_patch, nt, I_init, alpha, strain2_delay)

  if (!is.null(seed)) set.seed(seed)

  makefun <- get(sprintf("make_simulator_%s", platform))
  runfun <- get(sprintf("run_simulator_%s", platform))
  convfun <- get(sprintf("conv_%s", platform))

  ## cache the built simulator per-worker: first call builds, subsequent calls reuse.
  ## <<- assigns into FUN's closure env (discrete_run's frame locally; worker's
  ## deserialized closure copy in parallel), giving per-worker caching for free.
  mod_cache <- NULL

  FUN <- function(i) {
    if (is.null(mod_cache)) mod_cache <<- do.call(makefun, args)
    convfun(runfun(mod_cache))
  }

  if (nsim == 1) return(FUN(1))

  created_cl <- is.null(cl)
  cl <- cl %||% makeCluster(ncores)
  if (created_cl) {
    on.exit(stopCluster(cl))
    clusterSetRNGStream(cl)
  }

  clusterExport(cl, varlist = "FUN", envir = environment())
  clusterEvalQ(cl, { library(odin); library(dde); library(macpan2) })
  return(parLapply(cl = cl, X = seq.int(nsim), fun = FUN))

}

## discrete_run(nsim = 10)

sim_to_sum <- function(x, return_type = c("long", "wide")) {
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
                  type = sub("_.*", "", state))
}

