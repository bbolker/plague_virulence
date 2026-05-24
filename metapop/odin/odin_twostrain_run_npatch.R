library(odin)
library(dde) ## odin insists on this
library(tidyr)
library(dplyr)
library(ggplot2)
library(parallel)
library(patchwork)

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
                          n_strain = 2,
                          nt = 1000,
                          alpha = 1e-3,
                          strain2_delay = 0,
                          I_init = 10,
                          seed = NULL,
                          nsim = 1,
                          cl = NULL,
                          ncores = 1) {

  ## patch-level parameters (vectors, length n_patch)
  K_vec <- rep(K, length.out = n_patch)
  r_vec <- rep(r, length.out = n_patch)
  I_init <- rep(I_init, length.out = n_strain)
  if (!is.null(seed)) set.seed(seed)
  
  if (n_strain != 2) {
    stop("code not set up for n_strain > 2 (and not tested/may not work for n_strain < 2")
  }

  ## strain x patch parameters: matrices [n_patch, n_strain], rows = patches, cols = strains
  I_ini_mat <- matrix(rpois(n_strain*n_patch, lambda = I_init),
                      byrow = TRUE,
                      ncol = n_strain)

  S_ini_vec <- K_vec - rowSums(I_ini_mat)

  odin_file <- here::here("metapop/odin", "odin_twostrain_npatch.R")

  args <- tibble::lst(beta    = beta_vec,
                      r       = r_vec,
                      K       = K_vec,
                      I_ini   = I_ini_mat,
                      S_ini   = S_ini_vec,
                      alpha,
                      strain2_delay,
                      n_patch)

  FUN <- function(i) {
    gen_local <- suppressMessages(odin::odin(odin_file))
    mod_local <- do.call(gen_local$new, args)
    tt <- system.time(res <- mod_local$run(1:nt))
    attr(res, "time") <- tt
    res
  }

  if (nsim == 1) return(FUN(1))

  cl <- cl %||% { on.exit(stopCluster(cl)); cl <- makeCluster(ncores) }

  clusterExport(cl, varlist = c("odin_file", "args", "nt", "FUN"), envir = environment())
  clusterEvalQ(cl, { library(odin); library(dde) })
  clusterSetRNGStream(cl)
  return(parLapply(cl = cl, X = seq.int(nsim), fun = FUN))

}

sim_to_long <- function(x) {
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

sim_to_sum <- function(x, return_type = c("long", "wide")) {
  return_type <- match.arg(return_type)
  ## total population sizes over time
  total_pops <- rowSums(x)
  i1 <- x[, grepl(",1\\]", colnames(x))]
  i2 <- x[, grepl(",2\\]", colnames(x))]
  infected_patches <- cbind(I1 = rowSums(i1 > 0),
                       I2 = rowSums(i2 > 0))
  total_inf <- cbind(rowSums(i1), rowSums(i2))
  ret <- data.frame(total_pops, infected_patches, total_inf) |>
    as_tibble() |>
    setNames(c("S_pop", "I1_patch", "I2_patch", "I1_pop", "I2_pop"))
  if (return_type == "wide") return(ret)
  ret <- ret |>
    mutate(step = seq.int(nrow(ret))) |> 
    tidyr::pivot_longer(cols = -step, names_to = "state") |>
    mutate(var = sub(".*_", "", state),
           type = sub("_.*", "", state))
  ret
}

run1 <- run_twostrain(seed = 101)
cc1 <- sim_to_long(run1)
cc2 <- sim_to_sum(run1)

run_sep <- run_twostrain(seed = 101, alpha = 0)
cc_sep <- sim_to_long(run_sep)

run_twostrain(seed = 101, nsim = 20)

gg1 <- ggplot(cc1, aes(step, value, colour = state)) +
  geom_line(aes(group = interaction(state, patch)))  +
  scale_y_log10()

print(gg1)

print(gg_sep <- gg1 + cc_sep)

run_1strain <- run_twostrain(seed = 101, I_init = c(10, 0))
cc_1strain <- sim_to_long(run_1strain)
print(gg_1strain <- gg1 + cc_1strain)

ggsave(filename = "odin_twostrain_run_patch.png", plot = gg1)

ggplot(cc2, aes(step, value, colour = type)) +
  geom_line() +
  facet_wrap(~var, scale = "free") +
  scale_y_log10()
