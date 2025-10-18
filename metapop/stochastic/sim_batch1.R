source("simulation_funs.R")
dd <- expand.grid(R0 = seq(1.5, 3, by = 0.1),
                  alphavec = 5*10^seq(-6,-3, by = 0.5))

fn <- "sim_batch1_raw.rds"
res <- list()
cat("n R0 alpha total_pops inf_patch total_inf\n")
for (i in 1:nrow(dd)) {
  params <- params0
  params[["R0"]] <- dd$R0vec[i]
  params[["alpha"]] <- dd$alphavec[i]
  res[[i]] <- mult_sim_mp(nsim = 120, ncores = 12, params = params, seed = 100 + i)
  cat(i, unlist(dd[i,]), res[[i]], "\n")
  saveRDS(res, fn)
}

res <- lapply(res, sumfun1)
dd <- data.frame(dd, do.call(rbind, res))
saveRDS(dd, "sim_batch1.rds")
