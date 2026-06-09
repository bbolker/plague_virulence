library(plagueMetapop)
library(future)
library(dplyr)
library(ggplot2); theme_set(theme_bw())

## lower right panel of PIP plot, purple area (resident fine, invader loses)
R01 <- 2.5
R02 <- 1.5
K <- 1e5
alpha <- 1e-3
n_patch   <- 200
n_sim     <- 25L
dt        <- 0.1
stop_cond <- stop_both_extinct
strain2_delay <- round(100 / dt)
nt <- round(500/dt)
I_init <- cbind(rep(10, n_patch),
                c(10, rep(0, n_patch - 1L)))
runs <- discrete_run(beta_vec      = c(R01, R02),
                     K             = K,
                     r             = 0.125,
                     n_patch       = n_patch,
                     nt            = round(500 / dt),
                     alpha         = alpha,
                     I_init        = I_init,
                     gamma         = c(1, 1),
                     dt            = dt,
                     def_file      = "euler_odin_def.R",
                     strain2_delay = strain2_delay,
                     stop_cond     = stop_cond,
                     nsim          = n_sim,
                     platform      = "odin")

runs_x <- sum_runs(runs) |> dplyr::filter(var == "pop") |>
  mutate(across(value, ~ . /n_patch))

ggplot(runs_x, aes(step, value, colour = type)) +
  geom_line(aes(group = interaction(run, type)), alpha = 0.4) +
  scale_y_log10()

sumfun_discrete(runs, strain2_delay = strain2_delay)

