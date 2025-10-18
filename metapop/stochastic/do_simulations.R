source("simulation_funs.R")
system.time(res <- mult_sim_mp(nsim = 10, ncores = 10, params = params0, verbose = TRUE))
plotfun1(res)

dd <- expand.grid(R0 = seq(1.5, 3, by = 0.1),
                  alphavec = 5*10^seq(-6,-3, by = 0.5))

fn <- "sim_batch1.rds"
res <- list()
cat("n R0 alpha total_pops inf_patch total_inf\n")
for (i in 1:nrow(dd)) {
  params <- params0
  params[["R0"]] <- dd$R0vec[i]
  params[["alpha"]] <- dd$alphavec[i]
  res0 <- mult_sim_mp(nsim = 80, ncores = 12, params = params, seed = 100 + i)
  res[[i]] <- sumfun1(res0)
  cat(i, unlist(dd[i,]), res[[i]], "\n")
  saveRDS(res, fn)
}

dd <- data.frame(dd, do.call(rbind, res))
