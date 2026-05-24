library(tinytest)
using("tinysnapshot")

source("discrete_odin.R")
mod <- make_simulator_odin()
set.seed(101)
run <- run_simulator_odin(mod)
out <- conv_odin(run) |> dplyr::arrange(step, patch, state)

gg1 <- ggplot(out, aes(step, value, colour = state)) +
  geom_line(aes(group = interaction(state, patch)))  +
  scale_y_log10()


expect_snapshot_print(out, label = "odin_run1")

source("discrete_macpan2.R")
mod <- make_simulator_macpan2()
set.seed(101)
run <- run_simulator_macpan2(mod)
out <- conv_macpan2(run) |> dplyr::arrange(step, patch, state)

expect_snapshot_print(out, label = "macpan2_run1")

source("discrete_pureR.R")
set.seed(101)
pureR_run1 <- run_simulator_pureR()

expect_snapshot_print(out, label = "pureR_run1")
