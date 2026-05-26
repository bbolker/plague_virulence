library(plagueMetapop)
library(future)
library(ggplot2); theme_set(theme_bw())

nsim <- 100L

plan(multicore(workers = 10L))
set.seed(101)

run <- discrete_run(beta_vec = c(4, 0),
                    K        = 1e6,
                    r        = 0.125,
                    n_patch  = 1,
                    nt       = 500,
                    alpha    = 0,
                    I_init    = c(10, 0),
                    stop_cond = NULL,
                    nsim      = nsim,
                    platform  = "odin")

runx <- bind_rows(run, .id = "run")

dplyr::filter(runx, state == "I1") |>
  ggplot(aes(step, value)) +
  geom_line(aes(group = run)) +
  scale_y_log10()


