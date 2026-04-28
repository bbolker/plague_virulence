library("tinytest")
using("tinysnapshot")

source("sim_fun.R")
## debug(simulate_metapopulation_2strain)
expect_snapshot_print(simulate_metapopulation_2strain(n_patches = 2,
                                                      n_years  = 10,
                                                      R01 = 4, R02 = 4,
                                                      alpha = 3.2e-4,
                                                      rho = 8,
                                                      r = 0.5),
                      label = "sim1")
