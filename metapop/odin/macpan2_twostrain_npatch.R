library(macpan2)
library(ggplot2)
library(dplyr)

##' @param beta_vec length-2 vector of per-capita transmission rates
##' @param K carrying capacity (scalar or vector of length n_patch)
##' @param r host growth rate per disease generation (ditto)
##' @param n_patch number of patches
##' @param n_strain number of strains [only 2 supported; use beta_vec=c(b,0) for 1-strain]
##' @param nt time steps
##' @param alpha between-patch transmission probability
##' @param I_init mean initial infected per patch drawn from Poisson (length 1 or 2)
##' @param seed PRNG seed
##' @return long-format data frame with columns step, state, patch, value
run_twostrain_macpan2 <- function(
  beta_vec  = c(1.5, 2.5),
  K         = 1e4,
  r         = 0.125,
  n_patch   = 100,
  n_strain  = 2,
  nt        = 1000,
  alpha     = 1e-3,
  I_init    = 10,
  seed      = NULL
) {
  if (n_strain != 2) stop("only n_strain = 2 is supported")

  K_vec   <- matrix(rep(K, length.out = n_patch), ncol = 1)
  r_vec   <- matrix(rep(r, length.out = n_patch), ncol = 1)
  ones    <- matrix(1, nrow = n_patch, ncol = 1)
  I_init2 <- rep(I_init, length.out = 2)

  if (!is.null(seed)) set.seed(seed)
  I1_ini  <- matrix(rpois(n_patch, I_init2[1]), ncol = 1)
  I2_ini  <- matrix(rpois(n_patch, I_init2[2]), ncol = 1)
  S_ini   <- K_vec - I1_ini - I2_ini

  simple_sims(
    iteration_exprs = list(
      ## per-capita hazard of infection by each strain (mass action / K)
      hazard1     ~ beta1 * I1 / K,
      hazard2     ~ beta2 * I2 / K,
      sum_hazard  ~ hazard1 + hazard2,
      ## total infection probability (complement of survival probability)
      p_all       ~ 1 - exp(-sum_hazard),
      ## fraction of new infections that are strain 1; tiny avoids 0/0 when both absent
      p_SI1       ~ hazard1 / (sum_hazard + tiny),
      ## draw new infections then allocate between strains (sequential binomial)
      tot_incidence ~ rbinom(S, p_all),
      n_SI1       ~ rbinom(tot_incidence, p_SI1),
      n_SI2       ~ tot_incidence - n_SI1,
      ## logistic vital dynamics; clamp() prevents rpois() receiving negative argument
      ## when S > K (deaths in that regime are not modelled)
      delta_log   ~ r * S * (1 - S / K),
      pop_change  ~ rpois(clamp(delta_log)),
      ## colonization: strain-specific mean rate averaged over all patches,
      ## then n_patch independent Poisson draws (same rate for every patch)
      foi1        ~ (alpha / n_patch_s) * sum(I1),
      foi2        ~ (alpha / n_patch_s) * sum(I2),
      immig1      ~ rpois(foi1 * ones),
      immig2      ~ rpois(foi2 * ones),
      ## state updates
      S           ~ S  - tot_incidence + pop_change,
      I1          ~ n_SI1 + immig1,
      I2          ~ n_SI2 + immig2
    ),
    time_steps = nt,
    mats = list(
      ## state variables
      S             = S_ini,
      I1            = I1_ini,
      I2            = I2_ini,
      ## parameters
      beta1         = matrix(beta_vec[1]),
      beta2         = matrix(beta_vec[2]),
      K             = K_vec,
      r             = r_vec,
      alpha         = matrix(alpha),
      n_patch_s     = matrix(n_patch),
      ones          = ones,
      tiny          = matrix(1e-20),
      ## intermediates: must be in mats; all returned by simple_sims, filtered below
      hazard1       = matrix(0, n_patch, 1),
      hazard2       = matrix(0, n_patch, 1),
      sum_hazard    = matrix(0, n_patch, 1),
      p_all         = matrix(0, n_patch, 1),
      p_SI1         = matrix(0, n_patch, 1),
      tot_incidence = matrix(0, n_patch, 1),
      n_SI1         = matrix(0, n_patch, 1),
      n_SI2         = matrix(0, n_patch, 1),
      delta_log     = matrix(0, n_patch, 1),
      pop_change    = matrix(0, n_patch, 1),
      foi1          = matrix(0),
      foi2          = matrix(0),
      immig1        = matrix(0, n_patch, 1),
      immig2        = matrix(0, n_patch, 1)
    )
  ) |>
    filter(matrix %in% c("S", "I1", "I2")) |>
    rename(step = time, patch = row, state = matrix) |>
    select(step, state, patch, value)
}

run1 <- run_twostrain_macpan2(seed = 101)

gg1 <- ggplot(run1, aes(step, value, colour = state)) +
  geom_line(aes(group = interaction(state, patch))) +
  scale_y_log10() +
  theme_bw()

print(gg1)

ggsave(filename = "macpan2_twostrain_run_patch.png", plot = gg1)
