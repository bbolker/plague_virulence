library(tinytest)
using("tinysnapshot")

source(here::here("metapop/odin", "discrete_odin.R"))
mod <- make_simulator_odin()
set.seed(101)
run <- run_simulator_odin(mod)
out <- conv_odin(run) |> dplyr::arrange(step, patch, state)

run[nrow(run),]

library(ggplot2)
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

