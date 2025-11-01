source("simulation_funs.R")
retry <- TRUE
ncores <- 20
dd <- expand.grid(R0vec = seq(1.2, 3, by = 0.1),
                  alphavec = 5*10^seq(-6,-3, by = 0.5))

fn <- "sim_batch1_raw.rds"
res <- list()
start <- 1
if (retry && file.exists(fn)) {
  res <- readRDS(fn)
  start <- length(res) + 1
}
cat("n R0 alpha total_pops inf_patch total_inf\n")
for (i in start:nrow(dd)) {
  params <- params0
  params[["R0"]] <- dd$R0vec[i]
  params[["alpha"]] <- dd$alphavec[i]
  res[[i]] <- mult_sim_mp(nsim = 120, ncores = ncores, params = params, seed = 100 + i)
  cat(i, unlist(dd[i,]), sumfun1(res[[i]]), "\n")
  saveRDS(res, fn)
}

res <- lapply(res, sumfun1)
dd <- data.frame(dd, do.call(rbind, res))
saveRDS(dd, "sim_batch1.rds")
