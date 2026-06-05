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



plotfun <- function(run) {
    ii <- run |> dplyr::filter(state == "I1") |> dplyr::pull(value)    

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
}

