compile_odin <- function(def_file = "discrete_odin_def.R") {
  odin_file <- here::here("metapop/odin", def_file)
  gen <- suppressMessages(odin::odin(odin_file))
  attr(gen, "def_filename") <- def_file
  gen
}

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
  gen_local = NULL
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

  args <- tibble::lst(beta    = beta_vec,
                      r       = r_vec,
                      K       = K_vec,
                      I_ini   = I_ini_mat,
                      S_ini   = S_ini_vec,
                      I2_ini  = rep(0, n_patch),
                      alpha,
                      strain2_delay,
                      n_patch)

  if (grepl("^euler_", attr(gen_local, "def_filename"))) {
    args$gamma <- rep(gamma, length.out = n_strain)
    args$dt    <- dt
  }

  mod <- do.call(gen_local$new, args)
  attr(mod, "nt")        <- nt
  attr(mod, "gen")       <- gen_local
  attr(mod, "init_args") <- args
  mod
}

## @param stop_cond NULL (run all steps) or a function(row) -> logical called
##   on the last row of each chunk; return TRUE to stop early.
##   See stop_either_extinct() and stop_both_extinct() below.
## @param chunk number of steps per chunk when stop_cond is supplied
run_simulator_odin <- function(x, chunk = 50L, stop_cond = NULL) {
  nt        <- attr(x, "nt")
  gen       <- attr(x, "gen")
  init_args <- attr(x, "init_args")

  ## odin's $run() always restarts from initial conditions; it does NOT persist
  ## state between calls. drop_first removes the step-0 initial-conditions row.
  drop_first <- function(m) m[-1L, , drop = FALSE]

  if (is.null(stop_cond)) {
    return(drop_first(x$run(seq(0, nt))))
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
    out[[k]]  <- drop_first(res)
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

## Stopping conditions for stop_cond argument of run_simulator_odin().
## Operate on a raw odin output row (before conv_odin reshaping).
## Strain columns in the raw matrix are named "I[patch,strain]";
## grep on the trailing ",1]" / ",2]" selects each strain across all patches.
stop_either_extinct <- function(row) {
  sum(row[, grep(",1\\]$", colnames(row))]) == 0 ||
    sum(row[, grep(",2\\]$", colnames(row))]) == 0
}

stop_both_extinct <- function(row) {
  sum(row[, grep(",1\\]$", colnames(row))]) == 0 &&
    sum(row[, grep(",2\\]$", colnames(row))]) == 0
}

conv_odin <- function(x, format = c("long", "wide")) {
  format <- match.arg(format)
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

  

    
  
                                          
  
