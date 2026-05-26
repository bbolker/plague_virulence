library(tinytest)
using("tinysnapshot")

source(here::here("metapop/odin", "discrete_odin.R"))
mod <- make_simulator_odin()
set.seed(101)
run <- run_simulator_odin(mod)
out <- conv_odin(run) |> dplyr::arrange(step, patch, state)

run[nrow(run),]

library(ggplot2); theme_set(theme_bw())
gg1 <- ggplot(out, aes(step, value, colour = state)) +
  geom_line(aes(group = interaction(state, patch)))  +
  scale_y_log10()
print(gg1)

## expect_snapshot_print(out, label = "odin_run1")

source(here::here("metapop/odin", "discrete_macpan2.R"))
mod <- make_simulator_macpan2()
set.seed(101)
run <- run_simulator_macpan2(mod)
out <- conv_macpan2(run) |> dplyr::arrange(step, patch, state)

## expect_snapshot_print(out, label = "macpan2_run1")

afun <- function(x) dplyr::arrange(x, step, patch, state)

source(here::here("metapop/odin", "discrete_pureR.R"))
mod <- make_simulator_pureR()
set.seed(101)
run <- run_simulator_pureR(mod)
out <- conv_pureR(run)

## expect_snapshot_print(out, label = "pureR_run1")


source(here::here("metapop/odin", "discrete_run.R"))
run1 <- discrete_run(seed = 101)

run_sep <- discrete_run(seed = 101, alpha = 0)

gg1 <- ggplot(run1, aes(step, value, colour = state)) +
  geom_line(aes(group = interaction(state, patch)))  +
  scale_y_log10()

print(gg1)

print(gg_sep <- gg1 + run_sep)

run_1strain <- discrete_run(seed = 101, I_init = c(10, 0))
print(gg_1strain <- gg1 + run_1strain)


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

d1 <- d0 |> purrr::map_dfr(
  ~ . |> filter(state == "I1") |> select(-c(patch, state)),
  .id = "run")

ggplot(d1, aes(step, value)) + geom_line(aes(group = run)) +
  scale_x_log10() +
  scale_y_log10()
