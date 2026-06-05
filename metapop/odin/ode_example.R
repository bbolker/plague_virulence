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

ii <- run1 |> dplyr::filter(state == "I1") |> dplyr::pull(value)
plot(ii, log="y")

t1 <-traj_stats_ode(run1) |>
    cbind() |>
    as.data.frame() |>
    setNames("value") |>
    tibble::rownames_to_column("var")
    
ggplot(run1, aes(step, value)) + geom_line(aes(colour = state)) +
    scale_y_log10() +
    geom_hline(data = dplyr::filter(t1, !grepl("^t_", var)),
               aes(yintercept = value), lty = 2) +
    geom_vline(data = dplyr::filter(t1, grepl("^t_", var)),
               aes(xintercept = value), lty = 2)

cc <- attr(run1, "call")

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
