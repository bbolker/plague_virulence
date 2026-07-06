##' @rdname make_simulator_odin
##' @param n_strain number of strains [only 2 supported]
##' @export
make_simulator_pureR <- function(
  beta_vec  = c(1.5, 2.5),
  K         = 1e4,
  r         = 0.125,
  n_patch   = 100,
  n_strain  = 2,
  nt        = 1000,
  alpha     = 1e-3,
  I_ini_mat = make_I_ini_mat(10, n_patch),
  strain2_delay = 0
  ) {

  if (n_strain != 2) stop("only n_strain = 2 is supported")

  if (strain2_delay != 0) stop("strain2_delay not yet implemented for pureR")
  K_vec <- rep(K, length.out = n_patch)
  r_vec <- rep(r, length.out = n_patch)
  I     <- I_ini_mat
  S     <- K_vec - rowSums(I_ini_mat)

  ret <- list(run = function(nt) {

    ## pre-allocate output storage
    S_out <- matrix(NA_real_, nrow = nt, ncol = n_patch)
    I_out <- array(NA_real_,  dim = c(nt, n_patch, n_strain))

    for (t in seq_len(nt)) {
      ## hazard: n_patch x 2; sweep beta over columns, K over rows
      hazard     <- sweep(sweep(I, 2, beta_vec, `*`), 1, K_vec, `/`)
      sum_hazard <- rowSums(hazard)

      ## infection draws (sequential binomial, matching odin)
      p_all         <- 1 - exp(-sum_hazard)
      p_SI1         <- hazard[, 1] / (sum_hazard + 1e-20)
      tot_incidence <- rbinom(n_patch, S, p_all)
      n_SI1         <- rbinom(n_patch, tot_incidence, p_SI1)
      n_SI          <- cbind(n_SI1, tot_incidence - n_SI1)

      ## vital dynamics: rpois when S <= K, rbinom deaths when S > K
      delta_log  <- r_vec * S * (1 - S / K_vec)
      pos        <- delta_log >= 0
      pop_change <- integer(n_patch)
      if (any(pos))
        pop_change[ pos] <- rpois(sum(pos), delta_log[pos])
      if (any(!pos))
        pop_change[!pos] <- -rbinom(sum(!pos), S[!pos],
                                    pmin(1, -delta_log[!pos] / S[!pos]))

      ## colonization: n_patch independent Poisson draws per strain, same rate per strain
      foi   <- alpha * colSums(I) / n_patch
      immig <- matrix(rpois(n_patch * n_strain, rep(foi, each = n_patch)),
                      nrow = n_patch, ncol = n_strain)

      ## state updates
      S <- S - tot_incidence + pop_change
      I <- n_SI + immig

      S_out[t, ]   <- S
      I_out[t, , ] <- I
    }

    ## reshape to long format matching odin/macpan2 output
    dimnames(S_out) <- list(step = seq_len(nt), patch = seq_len(n_patch))
    dimnames(I_out) <- list(step = seq_len(nt), patch = seq_len(n_patch),
                            strain = seq_len(n_strain))

    S_long <- as.data.frame.table(S_out, responseName = "value")
    S_long$state <- "S"
    S_long$step  <- as.integer(S_long$step)
    S_long$patch <- as.integer(S_long$patch)

    I_long <- as.data.frame.table(I_out, responseName = "value")
    I_long$state <- paste0("I", I_long$strain)
    I_long$step  <- as.integer(I_long$step)
    I_long$patch <- as.integer(I_long$patch)

    ret <- rbind(S_long[, c("step", "state", "patch", "value")],
                 I_long[, c("step", "state", "patch", "value")])

    dplyr::as_tibble(ret)
  })
  environment(ret$run) <- environment()
  attr(ret, "nt") <- nt
  ret
}

##' @rdname make_simulator_odin
##' @export
run_simulator_pureR <- function(x) {
  nt <- attr(x, "nt")
  x$run(nt)
}

##' @rdname make_simulator_odin
##' @export
conv_pureR <- function(x) dplyr::as_tibble(x)
