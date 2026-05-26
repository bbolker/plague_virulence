##' Build a macpan2 simulator for the two-strain n-patch model.
##' Simulator creation is separated from trajectory generation so the two
##' steps can be timed independently.
##'
##' @inheritParams make_simulator_odin
##' @param seed PRNG seed (used for initial condition draws in R)
##' @return a macpan2 simulator object
##' @export
make_simulator_macpan2 <- function(
  beta_vec  = c(1.5, 2.5),
  K         = 1e4,
  r         = 0.125,
  n_patch   = 100,
  nt        = 1000,
  alpha     = 1e-3,
  I_init    = 10,
  strain2_delay = 0,
  seed      = NULL
) {

  n_strain <- 2 ## hard-coded

  if (strain2_delay != 0) stop("strain2_delay not yet implemented for macpan2")
  K_vec   <- matrix(rep(K, length.out = n_patch), ncol = 1)
  r_vec   <- matrix(rep(r, length.out = n_patch), ncol = 1)
  ones    <- matrix(1, nrow = n_patch, ncol = 1)
  I_init2 <- rep(I_init, length.out = 2)

  if (!is.null(seed)) set.seed(seed)
  I_ini <- matrix(rpois(n_strain * n_patch, lambda = I_init2),
                  nrow = n_patch, ncol = n_strain, byrow = TRUE)
  S_ini <- K_vec - rowSums(I_ini)

  spec <- macpan2::mp_tmb_model_spec(
    during = list(
      ## per-capita hazard: n_patch x 2; beta (1x2) and K (n_patch x 1) broadcast
      hazard        ~ I * beta / K,
      sum_hazard    ~ row_sums(hazard),
      ## total infection probability (complement of survival probability)
      p_all         ~ 1 - exp(-sum_hazard),
      ## fraction of new infections that are strain 1; tiny avoids 0/0 when both absent
      p_SI1         ~ row_sums(hazard * strain1_sel) / (sum_hazard + tiny),
      ## draw new infections then allocate between strains (sequential binomial)
      tot_incidence ~ rbinom(S, p_all),
      n_SI1         ~ rbinom(tot_incidence, p_SI1),
      ## build n_SI as n_patch x 2
      n_SI          ~ cbind(n_SI1, tot_incidence - n_SI1),
      ## logistic vital dynamics; clamp() prevents rpois() receiving negative argument
      ## when S > K (deaths in that regime are not modelled)
      delta_log     ~ r * S * (1 - S / K),
      pop_change    ~ rpois(clamp(delta_log)),
      ## colonization: foi is 1x2 (col_sums); ones %*% foi broadcasts to n_patch x 2
      foi           ~ (alpha / n_patch_s) * col_sums(I),
      immig         ~ rpois(ones %*% foi),
      ## state updates
      S             ~ S - tot_incidence + pop_change,
      I             ~ n_SI + immig
    ),
    default = list(
      ## state variables
      S             = S_ini,
      I             = I_ini,
      ## parameters
      beta          = matrix(beta_vec, nrow = 1),   # 1 x 2
      K             = K_vec,
      r             = r_vec,
      alpha         = matrix(alpha),
      n_patch_s     = matrix(n_patch),
      ones          = ones,
      tiny          = matrix(1e-20),
      strain1_sel   = matrix(c(1, 0), nrow = 1),   # 1 x 2: extracts strain 1 column
      ## intermediates (pre-allocated)
      hazard        = matrix(0, n_patch, 2),
      sum_hazard    = matrix(0, n_patch, 1),
      p_all         = matrix(0, n_patch, 1),
      p_SI1         = matrix(0, n_patch, 1),
      tot_incidence = matrix(0, n_patch, 1),
      n_SI1         = matrix(0, n_patch, 1),
      n_SI          = matrix(0, n_patch, 2),
      delta_log     = matrix(0, n_patch, 1),
      pop_change    = matrix(0, n_patch, 1),
      foi           = matrix(0, 1, 2),
      immig         = matrix(0, n_patch, 2)
    )
  )

  ## outputs = c("S","I") means only state variables are returned by mp_trajectory(),
  ## not all intermediates (unlike simple_sims which returns everything)
  macpan2::mp_simulator(spec, time_steps = nt, outputs = c("S", "I"))
}

##' Run a macpan2 simulator
##' @param x macpan2 simulator from make_simulator_macpan2()
##' @export
run_simulator_macpan2 <- function(x) {
  macpan2::mp_trajectory(x)
}

##' Reshape mp_trajectory() output to the long format used by odin/pureR versions
##' @param traj data frame returned by run_simulator_macpan2()
##' @export
conv_macpan2 <- function(traj) {
  traj |>
    dplyr::as_tibble() |>
    dplyr::filter(matrix %in% c("S", "I")) |>
    dplyr::rename(step = time, patch = row, state = matrix) |>
    ## account for zero-indexing of matrices
    dplyr::mutate(state = dplyr::if_else(state == "I", paste0("I", col + 1L), state)) |>
    dplyr::select(step, state, patch, value)
}
