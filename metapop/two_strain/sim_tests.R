library("tinytest")
using("tinysnapshot")

source("sim_fun.R")
debug(simulate_metapopulation_2strain)
sim1 <- simulate_metapopulation_2strain(n_patches = 2,
                                n_years  = 10,
                                R01 = 4, R02 = 4,
                                alpha = 3.2e-4,
                                rho = 8,
                                r = 0.5,
                                seed = 101)
expect_snapshot_print(sim1, label = "sim1")

load("polyfit.rda")
p1 <-pred_outcomes_poly(1.5, 2.0, 0.2, 0.1)
p2 <- pred_outcomes_poly(2.0, 1.5, 0.1, 0.2)
expect_equal(p1[1], p2[1])
expect_equal(p1[2], 1-p2[2])
expect_snapshot_print(p1, label = "pred1")
