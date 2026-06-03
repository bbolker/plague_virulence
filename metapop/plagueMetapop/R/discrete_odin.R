##' Compile an odin model definition file from inst/odin/
##' @param def_file filename of the odin DSL file in inst/odin/
##' @return an odin generator with a `def_filename` attribute
##' @export
compile_odin <- function(def_file = "discrete_odin_def.R") {
  odin_file <- system.file("odin", def_file, package = "plagueMetapop")
  gen <- suppressMessages(odin::odin(odin_file))
  attr(gen, "def_filename") <- def_file
  gen
}

##' Two-strain n-patch simulator interface (odin, macpan2, pureR backends)
##' @description
##' Three-function interface shared across all backends:
##' \itemize{
##'   \item \code{make_simulator_*()} constructs a simulator object
##'   \item \code{run_simulator_*()} runs it and returns raw output
##'   \item \code{conv_*()} converts raw output to a common long-format tibble
##'     with columns \code{step}, \code{state}, \code{patch}, \code{value}
##' }
##' @param beta_vec length-2 vector of per-capita transmission rates
##' @param K carrying capacity (scalar or vector of length n_patch)
##' @param r host growth rate per disease generation (ditto)
##' @param n_patch number of patches
##' @param nt time steps
##' @param alpha between-patch transmission probability
##' @param I_init mean initial infected per patch (Poisson draw; length 1 or 2)
##' @param strain2_delay steps before strain 2 is seeded (0 = immediate)
##' @param gamma length-2 vector of per-strain recovery rates (euler models only)
##' @param dt time-step size in disease-generation units (euler models only)
##' @param def_file odin DSL filename in inst/odin/
##' @param gen_local pre-compiled odin generator; compiled via compile_odin() if NULL
##' @return an odin model instance with nt, gen, and init_args attributes
##' @export
make_simulator_odin <- function(
  beta_vec  = c(1.5, 2.5),
  K         = 1e4,
  r         = 0.125,
  n_patch   = 100,
  nt        = 1000,
  alpha     = 1e-3,
  I_init    = 10,
  strain2_delay = 0,
  gamma     = c(1, 1),
  dt        = 1,
  def_file  = "discrete_odin_def.R",
  gen_local = NULL,
  logistic_growth = 1,
  reedfrost = 0
  ) {

  n_strain <- 2 ## hard-coded on purpose
  gen_local <- gen_local %||% compile_odin(def_file)

  ## patch-level parameters (vectors, length n_patch)
  K_vec <- rep(K, length.out = n_patch)
  r_vec <- rep(r, length.out = n_patch)
  I_init <- rep(I_init, length.out = n_strain)

  ## strain x patch parameters: matrices [n_patch, n_strain], rows = patches, cols = strains
  ## single-patch: use I_init directly (deterministic); multi-patch: Poisson draws
  I_ini_mat <- if (n_patch == 1L) {
    matrix(round(I_init), byrow = TRUE, ncol = n_strain)
  } else {
    matrix(rpois(n_strain * n_patch, lambda = I_init), byrow = TRUE, ncol = n_strain)
  }

  S_ini_vec <- K_vec - rowSums(I_ini_mat)

  def_fn <- attr(gen_local, "def_filename")

  ## Base args shared by all models
  args <- tibble::lst(beta  = beta_vec,
                      r     = r_vec,
                      K     = K_vec,
                      I_ini = I_ini_mat,
                      S_ini = S_ini_vec,
                      alpha,
                      n_patch)

  if (def_fn == "ode_odin_def.R") {
    ## ODE model: no step-based seeding; rates per unit time; no dt or reedfrost
    args <- c(args, tibble::lst(gamma = rep(gamma, length.out = n_strain),
                                logistic_growth))
  } else {
    ## Discrete and Euler models: step-based strain-2 seeding via I2_ini/strain2_delay
    args$I2_ini        <- rep(0, n_patch)
    args$strain2_delay <- strain2_delay
    if (grepl("^euler_", def_fn)) {
      args <- c(args, tibble::lst(gamma = rep(gamma, length.out = n_strain), dt))
      if (def_fn == "euler_odin_def.R")
        args <- c(args, tibble::lst(logistic_growth, reedfrost))
    }
  }

  mod <- do.call(gen_local$new, args)
  attr(mod, "nt")        <- nt
  attr(mod, "gen")       <- gen_local
  attr(mod, "init_args") <- args
  mod
}

##' @rdname make_simulator_odin
##' @param x simulator object from the corresponding \code{make_simulator_*()} function
##' @param chunk steps per chunk; only used when stop_cond is non-NULL
##' @param stop_cond NULL (run all nt steps) or a function(row) -> logical
##'   called on the last row of each chunk; return TRUE to stop early.
##'   See stop_either_extinct() and stop_both_extinct().
##' @export
run_simulator_odin <- function(x, chunk = 50L, stop_cond = NULL) {
  nt        <- attr(x, "nt")
  gen       <- attr(x, "gen")
  init_args <- attr(x, "init_args")

  ## odin's $run() always restarts from initial conditions; it does NOT persist
  ## state between calls. drop_first removes the step-0 initial-conditions row.
  drop_first <- function(m) m[-1L, , drop = FALSE]

  if (is.null(stop_cond)) {
    return(x$run(seq(0, nt)))
  }

  ## Chunked early-stopping (approach 2):
  ## Each chunk runs seq(0, chunk_len) on a fresh model instance initialised
  ## from the previous chunk's last-row state. This works around odin resetting
  ## to initial conditions on every $run() call.
  ##
  ## strain2_delay is adjusted each chunk: once the seeding step has passed
  ## (orig_delay < steps_elapsed) we set it to .Machine$integer.max so it
  ## never fires again; otherwise we decrement by the elapsed step count so
  ## the seeding fires at the right absolute time.
  breaks             <- unique(c(seq(0L, nt, by = chunk), nt))
  out                <- vector("list", length(breaks) - 1L)
  mod                <- x
  orig_strain2_delay <- init_args$strain2_delay

  for (k in seq_along(out)) {
    chunk_len <- breaks[k + 1L] - breaks[k]
    res       <- mod$run(seq(0L, chunk_len))
    out[[k]]  <- if (k==1) res else drop_first(res)
    out[[k]][, 1L] <- out[[k]][, 1L] + breaks[k]   ## shift to absolute step numbers

    if (stop_cond(res[nrow(res), , drop = FALSE])) break

    if (k < length(out)) {
      t_elapsed <- breaks[k + 1L]
      last_row  <- res[nrow(res), , drop = FALSE]
      cn        <- colnames(last_row)

      new_args           <- init_args
      new_args$S_ini       <- last_row[, grepl("^S\\[",  cn)]
      new_args$I_ini[, 1L] <- last_row[, grepl(",1\\]$", cn)]
      new_args$I2_ini      <- last_row[, grepl(",2\\]$", cn)]
      new_args$strain2_delay <- if (orig_strain2_delay < t_elapsed) {
        .Machine$integer.max   ## seeding already fired; suppress forever
      } else {
        orig_strain2_delay - t_elapsed
      }

      mod <- do.call(gen$new, new_args)
    }
  }

  do.call(rbind, out[!vapply(out, is.null, logical(1L))])
}

##' Stopping conditions for run_simulator_odin()
##' @description
##'   \code{stop_either_extinct()} is a factory: call it once to create a
##'   stateful closure, then pass the closure as \code{stop_cond}.  The closure
##'   returns TRUE when strain 1 is globally extinct, or when strain 2 was
##'   previously present and is now globally extinct.  Strain 2 is considered
##'   absent (not yet introduced) until the closure first observes I2 > 0, so
##'   simulations with \code{strain2_delay > 0} are not stopped prematurely.
##'   Because \code{furrr} serialises the closure independently for each task,
##'   each simulation gets its own fresh \code{strain2_seen = FALSE}.
##'
##'   \code{stop_both_extinct(row)} is a plain function (not a factory) that
##'   returns TRUE only when both strains are simultaneously globally extinct.
##'   It does not need factory wrapping because I1 > 0 before strain-1 goes
##'   extinct, so the condition cannot fire before strain 2 is seeded.
##' @return \code{stop_either_extinct()} returns a \code{function(row)} closure.
##' @examples
##' ## pass the closure, not the factory:
##' ## runs <- discrete_run(nsim = 20, stop_cond = stop_either_extinct())
##' @export
stop_either_extinct <- function() {
  strain2_seen <- FALSE
  function(row) {
    I1 <- sum(row[, grep(",1\\]$", colnames(row))])
    I2 <- sum(row[, grep(",2\\]$", colnames(row))])
    if (I2 > 0) strain2_seen <<- TRUE
    I1 == 0 || (strain2_seen && I2 == 0)
  }
}

##' @rdname stop_either_extinct
##' @param row single-row matrix from odin$run(), column names "I[patch,strain]"
##' @export
stop_both_extinct <- function(row) {
  sum(row[, grep(",1\\]$", colnames(row))]) == 0 &&
    sum(row[, grep(",2\\]$", colnames(row))]) == 0
}

##' @rdname make_simulator_odin
##' @param format "long" (default) or "wide"
##' @export
conv_odin <- function(x, format = c("long", "wide")) {
  format <- match.arg(format)
  ## ODE models label the time column 't'; normalise to 'step' for uniform output
  colnames(x)[colnames(x) == "t"] <- "step"
  ## convert S[i] -> S_i, I[i,j] -> Ij_i
  cn <- colnames(x)
  cn <- cn |>
    gsub(pattern = "\\[([0-9]+),([0-9]+)", replacement = "\\2[\\1") |>
    gsub(pattern = "\\[([0-9]+)\\]", replacement = "_\\1")
  colnames(x) <- cn
  if (format == "wide") return(x)
  x <- x |>
    dplyr::as_tibble() |>
    tidyr::pivot_longer(-step, names_to = "state") |>
    dplyr::mutate(patch = as.numeric(sub(".*_", "", state)),
                  state = sub("_.*", "", state)) |>
    dplyr::select(step, state, patch, value)
  x
}
