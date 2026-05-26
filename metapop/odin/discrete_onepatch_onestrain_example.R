library(dplyr)
library(tidyr)
library(odin)
library(dde)
library(future)
library(ggplot2); theme_set(theme_bw())

source(here::here("metapop/odin", "discrete_run.R"))

nsim <- 100L

plan(multicore(workers = 10L))
set.seed(101)

run <- discrete_run(beta_vec = c(4, 0),
                    K        = 1e6,
                    r        = 0.125,
                    n_patch  = 1,
                    nt       = 500,
                    alpha    = 0,
                    I_init   = c(10, 0),
                    nsim     = nsim,
                    platform = "odin")

runx <- bind_rows(run, .id = "run")

dplyr::filter(runx, state == "I1") |>
  ggplot(aes(step, value)) +
  geom_line(aes(group = run)) +
  scale_y_log10()


  bind_cols(dd[i, ], as.data.frame(as.list(sumfun_discrete(runs))))
})
plan(sequential)

out <- bind_rows(result_list)
saveRDS(out, "discrete_onepatch_onestrain_extinct.rds")
