library(plagueMetapop)
library(future)
library(ggplot2); theme_set(theme_bw())

nsim <- 100L
dt   <- 0.2

plan(multicore(workers = 10L))
set.seed(101)

runs <- discrete_run(beta_vec  = c(4, 0),
                    K         = 1e5,
                    r         = 0.125,
                    n_patch   = 1,
                    nt        = round(100 / dt),
                    alpha     = 0,
                    I_init    = c(2, 0),
                    gamma     = 1,
                    dt        = dt,
                    def_file  = "euler_odin_def.R",
                    stop_cond = NULL,
                    nsim      = nsim,
                    platform  = "odin")
plan(sequential)

runx <- bind_rows(runs, .id = "run")

runx2 <- dplyr::filter(runx, state %in% c("I1", "S"))

ggplot(runx2, aes(step, value)) +
  geom_line(aes(group = interaction(run, state), colour = state))

ggsave("euler_onepatch_onestrain_example.pdf")
sumfun_discrete(runs)
