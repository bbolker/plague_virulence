library(plagueMetapop)
library(future)
library(furrr)

plan(multicore(workers = 15))

dt <- 0.01
run1 <- discrete_run(beta_vec  = c(4, 0),
                    K         = 1,
                    r         = 0.125,
                    n_patch   = 1,
                    nt        = round(10/ dt),
                    alpha     = 0,
                    I_init    = c(1e-5, 0),
                    I_ini_method = "fixed",
                    gamma     = 1,
                    dt        = dt,
                    def_file  = "ode_odin_def.R",
                    stop_cond = NULL,
                    nsim      = 1,
                    platform  = "odin") |>
    dplyr::filter(state != "I2")

cc <- attr(run1, "call")

dt <- 0.01
ng <- 51
dd <- expand.grid(beta = seq(1.1, 10, length = ng),
                  r = 10^seq(-3, log(0.5), length = ng))
res <- vector("list", length = 2) |> setNames(c("linear", "logistic"))
for (nm in names(res)) {
    logistic <- as.numeric(nm == "logistic")
    cc$logistic <- logistic
    cc_nm <- cc
    rows <- future_map(seq_len(nrow(dd)), function(i) {
        cc_i <- cc_nm
        cc_i$beta_vec <- c(dd$beta[i], 0)
        cc_i$r <- dd$r[i]
        ## need to re-evaluate here; lazy eval not quite smart enough
        cc_i$nt <- round(10 / dd$r[i] / dt)
        traj_stats_ode(eval(cc_i))
    }, .options = furrr_options(seed = TRUE))
    res[[nm]] <- dplyr::bind_cols(dd, dplyr::bind_rows(rows))
} ## loop over demography


res <- dplyr::bind_rows(res, .id = "demography")
saveRDS(res, "ode_trough.rds")
