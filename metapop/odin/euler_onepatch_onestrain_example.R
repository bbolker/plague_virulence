library(dplyr)
library(tidyr)
library(odin)
library(dde)
library(future)
library(ggplot2); theme_set(theme_bw())

source(here::here("metapop/odin", "discrete_run.R"))

nsim <- 1L
dt   <- 0.1

plan(multicore(workers = 10L))
set.seed(101)

run <- discrete_run(beta_vec  = c(2, 0),
                    K         = 1e4,
                    r         = 0.125,
                    n_patch   = 1,
                    nt        = round(100 / dt),
                    alpha     = 0,
                    I_init    = c(2, 0),
                    gamma     = 1,
                    dt        = dt,
                    def_file  = "euler_det_odin_def.R",
                    stop_cond = NULL,
                    nsim      = nsim,
                    platform  = "odin")
plan(sequential)

runx <- bind_rows(run, .id = "run")

runx2 <- dplyr::filter(runx, state %in% c("I1", "S"))

ggplot(runx2, aes(step, value)) +
  geom_line(aes(group = interaction(run, state), colour = state))
