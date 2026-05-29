library(plagueMetapop)
library(future)
library(ggplot2); theme_set(theme_bw())

nsim <- 20L
dt   <- 0.1

plan(multicore(workers = 20L))
set.seed(101)

runs <- discrete_run(beta_vec  = c(2, 2.5),
                    K         = 1e4,
                    r         = 0.125,
                    n_patch   = 200,
                    nt        = round(100 / dt),
                    alpha     = 10^-5.5,
                    I_init    = c(10, 10),
                    gamma     = 1,
                    strain2_delay = round(50/dt),
                    dt        = dt,
                    def_file  = "euler_odin_def.R",
                    stop_cond = stop_both_extinct,
                    nsim      = nsim,
                    platform  = "odin")
plan(sequential)

runx <- dplyr::bind_rows(runs, .id = "run")

runx2 <- dplyr::filter(runx, state %in% c("I1", "I2"),
                       step %% 10 == 1)

ggplot(runx2, aes(step, value)) +
  geom_line(aes(group = interaction(run, state, patch), colour = state),
            alpha = 0.4) +
  scale_colour_brewer(palette = "Dark2") +
  scale_y_log10()

ggsave("euler_twostrain_example.pdf")
sumfun_discrete(runs)
