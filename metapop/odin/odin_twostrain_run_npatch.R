library(odin)
library(dde) ## odin insists on this
library(tidyr)
library(dplyr)
library(ggplot2)

##' @param K carrying capacity (scalar or vector of length n_patch)
##' @param r host growth rate per disease generation (ditto)
##' @param n_patch number of patches
##' @param n_strains number of strains [HARD-CODED]
##' @param nt time steps
##' @param alpha between-patch transmission  probability
##' @param I_init initial infection (chosen as rpois)
##' @param seed PRNG seed
run_twostrain <- function(beta_vec = c(1.5, 2.5),
                          K = 1e4,
                          r = 0.125,
                          n_patch = 100,
                          n_strains = 2,
                          nt = 1000,
                          alpha = 1e-3,
                          I_init = 10,
                          seed = NULL) {

  ## patch-level parameters (vectors, length n_patch)
  K_vec <- rep(K, length.out = n_patch)
  r_vec <- rep(r, length.out = n_patch)

  if (!is.null(seed)) set.seed(seed)
  
  if (n_strains != 2) {
    stop("code not set up for n_strains > 2 (and not tested/may not work for n_strains < 2")
  }
  ## strain x patch parameters: matrices [n_patch, n_strains], rows = patches, cols = strains

  I_ini_mat <- matrix(rpois(n_strains*n_patch, lambda = I_init),
                      ncol = n_strains)

  S_ini_vec <- K_vec - rowSums(I_ini_mat)

  odin_file <- here::here("metapop/odin", "odin_twostrain_npatch.R")
  gen <- odin::odin(odin_file)

  args <- tibble::lst(beta    = beta_vec,
                      r       = r_vec,
                      K       = K_vec,
                      I_ini   = I_ini_mat,
                      S_ini   = S_ini_vec,
                      alpha,
                      n_patch)

  mod <- do.call(gen$new, args)

  tt <- system.time(res <- mod$run(1:nt))

  attr(res, "time") <- tt

  return(res)
}

conv_fun <- function(x) {
  ret <- as.data.frame.table(x[,!colnames(x) %in% "step"]) |>
    as_tibble() |>
    transmute(step = as.numeric(Var1), ## grab factor levels
              state = case_when(
                grepl("^S", Var2) ~ gsub("\\[", ".", gsub("\\]", "", Var2)),
                grepl("^I", Var2) ~ gsub("\\[([0-9]+),([0-9]+)\\]", "\\2.\\1", Var2)),
              value = Freq) |>
    mutate(patch = as.numeric(sub(".*\\.", "", state)),
           state = sub("\\..*", "", state)) |>
    arrange(step, patch)
  return(ret)
}

run1 <- run_twostrain(seed = 101)
cc1 <- conv_fun(run1)

run_sep <- run_twostrain(seed = 101, alpha = 0)
cc_sep <- conv_fun(run_sep)

gg1 <- ggplot(cc1, aes(step, value, colour = state)) +
  geom_line(aes(group = interaction(state, patch)))  +
  scale_y_log10()

print(gg1)

print(gg_sep <- gg1 + cc_sep)

ggsave(filename = "odin_twostrain_run_patch.png", plot = gg1)
