source("simulation_funs.R")
system.time(res <- mult_sim_mp(nsim = 5, params = params0, verbose = TRUE))

