make_simulator_odin <- function(
  beta_vec  = c(1.5, 2.5),
  K         = 1e4,
  r         = 0.125,
  n_patch   = 100,
  nt        = 1000,
  alpha     = 1e-3,
  I_init    = 10,
  strain2_delay = 0
  ) {

  n_strain <- 2 ## hard-coded on purpose
  odin_file <- here::here("metapop/odin", "discrete_odin_def.R")

  ## patch-level parameters (vectors, length n_patch)
  K_vec <- rep(K, length.out = n_patch)
  r_vec <- rep(r, length.out = n_patch)
  I_init <- rep(I_init, length.out = n_strain)
  
  ## strain x patch parameters: matrices [n_patch, n_strain], rows = patches, cols = strains
  I_ini_mat <- matrix(rpois(n_strain*n_patch, lambda = I_init),
                      byrow = TRUE,
                      ncol = n_strain)

  S_ini_vec <- K_vec - rowSums(I_ini_mat)

  args <- tibble::lst(beta    = beta_vec,
                      r       = r_vec,
                      K       = K_vec,
                      I_ini   = I_ini_mat,
                      S_ini   = S_ini_vec,
                      alpha,
                      strain2_delay,
                      n_patch)

  gen_local <- suppressMessages(odin::odin(odin_file))
  mod <- do.call(gen_local$new, args)
  attr(mod, "nt") <- nt
  mod
}

## @param stop_cond NULL (run all steps) or a function(row) -> logical called
##   on the last row of each chunk; return TRUE to stop early.
##   See stop_either_extinct() and stop_both_extinct() below.
## @param chunk number of steps per chunk when stop_cond is supplied
run_simulator_odin <- function(x, chunk = 50L, stop_cond = NULL) {
  nt <- attr(x, "nt")

  ## odin returns state at every requested step, including the first one in the
  ## vector. For consecutive chunks the shared boundary step appears as both the
  ## last row of one chunk and the first row of the next, so it would be
  ## duplicated in the combined output. Dropping the first row of every chunk
  ## with [-1, ] removes the duplicate; for the first chunk this also drops the
  ## step-0 initial-conditions row, matching the original single-call behaviour.
  drop_first <- function(m) m[-1L, , drop = FALSE]

  if (is.null(stop_cond)) {
    return(drop_first(x$run(seq(0, nt))))
  }

  breaks <- unique(c(seq(0L, nt, by = chunk), nt))
  out <- vector("list", length(breaks) - 1L)
  for (k in seq_along(out)) {
    res      <- x$run(seq(breaks[k], breaks[k + 1L]))
    out[[k]] <- drop_first(res)
    if (stop_cond(res[nrow(res), , drop = FALSE])) break
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

  

    
  
                                          
  
