library(plagueMetapop)
library(ggplot2); theme_set(theme_bw())

dt <- 0.01
run1 <- discrete_run(beta_vec  = c(4, 0),
                    K         = 1,
                    r         = 0.125,
                    n_patch   = 1,
                    nt        = round(50 / dt),
                    alpha     = 0,
                    I_init    = c(0.001, 0),
                    I_ini_method = "fixed",
                    gamma     = 1,
                    dt        = dt,
                    def_file  = "ode_odin_def.R",
                    stop_cond = NULL,
                    nsim      = 1,
                    platform  = "odin") |>
  dplyr::filter(state != "I2")

cc <- attr(run1, "call")
dt <- 0.005
run2 <- eval(cc) |>   dplyr::filter(state != "I2")

ggplot(run1, aes(step, value, colour = state)) + geom_line() +
  scale_y_log10() +
  geom_line(data = run2, lty = 2)

t1 <-traj_stats_ode(run1)
t2 <-traj_stats_ode(run2)
all.equal(t1, t2, tolerance = 1e-3)

dt <- 0.01
dd <- expand.grid(beta = seq(1.1, 10, length = 51),
                  r = 10^seq(-3, log(0.5), length = 51))
res <- list()
pb <- txtProgressBar(max = nrow(dd), style = 3)
for (i in 1:nrow(dd)) {
  setTxtProgressBar(pb, i)
  cc$beta_vec <- c(dd$beta[i], 0)
  cc$r <- dd$r[i]
  runx <- eval(cc)
  res[[i]] <- traj_stats_ode(runx)
}
rm(pb)
