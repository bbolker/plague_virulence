library(plagueMetapop)
library(future)

## basic odin run
mod <- make_simulator_odin()
set.seed(101)
run <- run_simulator_odin(mod)
out <- conv_odin(run) |> dplyr::arrange(step, patch, state)

## basic macpan2 run
mod <- make_simulator_macpan2()
set.seed(101)
run <- run_simulator_macpan2(mod)
out <- conv_macpan2(run) |> dplyr::arrange(step, patch, state)

## basic pureR run
mod <- make_simulator_pureR()
set.seed(101)
run <- run_simulator_pureR(mod)
out <- conv_pureR(run)

## discrete_run: single runs
run1      <- discrete_run(seed = 101)
run_sep   <- discrete_run(seed = 101, alpha = 0)
run_1strain <- discrete_run(seed = 101, I_init = c(10, 0))

## NOTE: n_patch=1 with default K=1e4 and beta=c(2,1) — strain 1 (R0=2) is
## unlikely to go extinct within nt=1000 steps for most runs, but strain 2
## is never introduced (I_init[2]=0, alpha=0) so ext_prob.I2==1 is guaranteed.
## The ext_prob.I1 < 1 test would fail if all 6 runs went extinct, which is
## extremely unlikely at this K but could happen at much smaller K.
plan(multicore(workers = 2))
d0 <- discrete_run(beta_vec = c(2, 1), I_init = c(10, 0),
                   n_patch = 1, nsim = 6)
plan(sequential)

expect_equal(length(d0), 6L,
             info = "discrete_run returns one result per simulation")
expect_true(all(sapply(d0, is.data.frame)),
            info = "each parallel result is a data frame")
s0 <- sumfun_discrete(d0)
expect_equal(s0[["ext_prob.I2"]], 1,
             info = "strain 2 absent throughout (I_init[2]=0, alpha=0)")
expect_true(s0[["ext_prob.I1"]] < 1,
            info = "strain 1 (beta=2) persists in at least some runs")

d1 <- d0 |> purrr::map(\(x) dplyr::filter(x, state == "I1") |>
                          dplyr::select(-c(patch, state))) |>
  dplyr::bind_rows(.id = "run")
