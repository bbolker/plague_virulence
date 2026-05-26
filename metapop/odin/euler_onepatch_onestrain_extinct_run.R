library(dplyr)
library(tidyr)
library(odin)
library(dde)
library(future)

source(here::here("metapop/odin", "discrete_run.R"))

nsim <- 100L
dt   <- 0.1

dd <- expand.grid(R0 = seq(1.1, 5, by = 0.1),
                  K  = 10^(seq(3, 6, by = 0.5)))

plan(multicore(workers = 14L))
set.seed(101)

result_list <- lapply(seq_len(nrow(dd)), function(i) {
  cat(dd$R0[i], dd$K[i], "\n")
  runs <- discrete_run(beta_vec = c(dd$R0[i], 0),
                       K        = dd$K[i],
                       r        = 0.125,
                       n_patch  = 1,
                       nt       = round(200 / dt),
                       alpha    = 0,
                       I_init   = c(10, 0),
                       gamma    = c(1, 1),
                       dt       = dt,
                       def_file  = "euler_odin_def.R",
                       stop_cond = NULL,
                       nsim      = nsim,
                       platform  = "odin")
  bind_cols(dd[i, ], as.data.frame(as.list(sumfun_discrete(runs))))
})
plan(sequential)

out <- bind_rows(result_list)
saveRDS(out, "euler_onepatch_onestrain_extinct.rds")
