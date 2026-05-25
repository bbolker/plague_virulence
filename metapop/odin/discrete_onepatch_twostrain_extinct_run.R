library(dplyr)
library(tidyr)
library(odin)
library(dde) ## odin insists on this
library(future)

source(here::here("metapop/odin", "discrete_run.R"))

nsim <- 20

first_zero <- function(x) which(x == 0)[1]

## parameter note
## plague generation time ~ 10 - 20 days
## rat $r$ ~ 3 / year?  ~ 0.125?
ctr <- 0 ## checkpoint counter
sumfun <- function(beta1, beta2, K = 1e6, I_init = c(10, 10)) {
  runs <- discrete_run(beta_vec = c(beta1, beta2), K = K, r = 0.125,
                       n_patch = 1, nt = 1000, alpha = 0, I_init = I_init,
                       nsim = nsim, platform = "odin")
  if (!is.list(runs)) runs <- list(runs)
  ext <- sapply(runs, function(traj) {
    c(I1 = first_zero(traj$value[traj$state == "I1"]),
      I2 = first_zero(traj$value[traj$state == "I2"]))
  })
  rowMeans(ext, na.rm = TRUE)
}

R0vec <- seq(1, 5, by = 0.025)
dd <- expand.grid(R01 = R0vec, R02 = R0vec)

plan(multicore(workers = 10))
set.seed(101)
system.time(
  res <- apply(dd, 1, \(x) sumfun(x[1], x[2]))
)
plan(sequential)

res_df <- t(res) |> as.data.frame() |> bind_cols(dd)
dd2    <- res_df |> pivot_longer(cols = c("I1", "I2"))
dd2_i1 <- dd2 |> dplyr::filter(name == "I1")

fn <- "discrete_onepatch_twostrain_extinct.rds"
saveRDS(dd2_i1, file = fn)
