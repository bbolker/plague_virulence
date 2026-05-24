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

run_simulator_odin <- function(x) {
  nt <- attr(x, "nt")
  ## run from step 0 so step==0 fires during the first transition;
  ## drop the initial-conditions row (step 0) to match other platforms
  res <- x$run(seq(0, nt))
  res[-1, , drop = FALSE]
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

  

    
  
                                          
  
