library(plagueMetapop)
library(ggplot2); theme_set(theme_bw())
library(dplyr)

dt <- 0.01
r <- 0.003 ## was 0.125
run1 <- discrete_run(beta_vec  = c(7, 0), ## was (4, 0)
                    K         = 1,
                    r         = r, ## was 0.125
                    n_patch   = 1,
                    nt        = round(10 / r / dt),
                    alpha     = 0,
                    I_init    = c(0.001, 0),
                    I_ini_method = "fixed",
                    dt        = dt,
                    def_file  = "ode_odin_def.R",
                    logistic_growth = 1,
                    nsim      = 1,
                    stop_cond = NULL,
                    platform  = "odin") |>
    dplyr::filter(state != "I2")

plotfun <- function(run) {
    ii <- run |> filter(state == "I1") |> pull(value)    

    t0 <-traj_stats_ode(run)
    beg <- t0[["t_enter.boundary"]]
    end <- t0[["t_eqS"]]
    dd_trough <- run |> filter(between(step, beg, end))
    get_var <- function(v) {
      vn <- deparse(substitute(v))
      filter(dd_trough, state==vn) |> pull(value)
    }
    trough_S <- get_var(S)
    trough_step <- unique(dd_trough$step)
    t1 <- t0 |>
      cbind() |>
      as.data.frame() |>
      setNames("value") |>
      tibble::rownames_to_column("var") |>
      filter(var != "trough_area") |>
      mutate(var2 = case_when(
               grepl("S", var) ~ "S",
               grepl("I", var) ~ "I1",
               grepl("bound|trough", var) ~ "trough"))
    gg1 <- ggplot(run1, aes(step, value)) +
      geom_line(aes(colour = state, linetype = state)) +
      scale_y_log10() +
      geom_hline(data = filter(t1, !grepl("^t_", var)),
                 aes(yintercept = value, colour = var2, linetype  = var2), lty = 2) +
      geom_vline(data = filter(t1, grepl("^t_", var)),
                 aes(xintercept = value, colour = var2, linetype = var2), lty = 2) +
      annotate(geom = "polygon",
               x = c(trough_step, rev(trough_step)),
               y = c(trough_S, rep(t0[["eq_S"]], length(trough_S))),
               colour = NA, fill = "black", alpha = 0.3) +
      scale_colour_brewer(palette = "Dark2")

    return(gg1)
}

run1 |> filter(state=="I1") |> pull(value) |> min()

plotfun(run1)
ggsave(plotfun(run1), file = "ode_trough_example.pdf")
