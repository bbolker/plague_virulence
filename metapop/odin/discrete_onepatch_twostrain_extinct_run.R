library(dplyr)
library(tidyr)
library(odin)
library(dde) ## odin insists on this
library(future)
library(optparse)

source(here::here("metapop/odin", "discrete_run.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option(c("-s", "--stepR0"), type = "double",  default = 0.025,
              help = "step size of the R0 grid [default %default]"),
  make_option(c("-n", "--nsim"),   type = "integer", default = 20L,
              help = "simulations averaged per R0 pair [default %default]")
)))

nsim <- opt$nsim

first_zero <- function(x) which(x == 0)[1]

## parameter note
## plague generation time ~ 10 - 20 days
## rat $r$ ~ 3 / year?  ~ 0.125?
ctr             <- 0
checkpoint_freq <- 100
checkpoint_file <- "discrete_onepatch_twostrain_extinct_checkpoint.rds"

sumfun <- function(beta1, beta2, K = 1e6, I_init = c(10, 10)) {
  runs <- discrete_run(beta_vec = c(beta1, beta2), K = K, r = 0.125,
                       n_patch = 1, nt = 1000, alpha = 0, I_init = I_init,
                       stop_cond = NULL, nsim = nsim, platform = "odin")
  if (!is.list(runs)) runs <- list(runs)
  ext <- sapply(runs, function(traj) {
    c(I1 = first_zero(traj$value[traj$state == "I1"]),
      I2 = first_zero(traj$value[traj$state == "I2"]))
  })
  result <- rowMeans(ext, na.rm = TRUE)
  ctr <<- ctr + 1
  res[, ctr] <<- result
  if (ctr %% checkpoint_freq == 0) saveRDS(res, checkpoint_file)
  result
}

R0vec <- seq(1, 5, by = opt$stepR0)
dd    <- expand.grid(R01 = R0vec, R02 = R0vec)

## run once on first row to determine output dimensions; 2-row placeholder
## lets the call succeed before res is properly initialised
res      <- matrix(NA_real_, 2L, nrow(dd))
test_out <- sumfun(dd[1, 1], dd[1, 2])
res      <- matrix(NA_real_, nrow = length(test_out), ncol = nrow(dd),
                   dimnames = list(names(test_out), NULL))
ctr      <- 0  ## reset; test run does not count toward the grid sweep

plan(multicore(workers = 10))
set.seed(101)
system.time(
  apply(dd, 1, \(x) sumfun(x[1], x[2]))
)
plan(sequential)

res_df <- t(res) |> as.data.frame() |> bind_cols(dd)
dd2    <- res_df |> pivot_longer(cols = c("I1", "I2"))
dd2_i1 <- dd2 |> dplyr::filter(name == "I1")

fn <- "discrete_onepatch_twostrain_extinct.rds"
saveRDS(dd2_i1, file = fn)
