source("simulation_funs.R")

fn_raw <- "outputs/sim_raw.rds"
fn_sum <- "outputs/sim.rds"

retry <- TRUE
ncores <- 16
nsim <- 120

dd <- expand.grid(
  R0vec = seq(1.2, 3.0, by = 0.2),
  kappavec = 10^seq(-6, -3, by = 0.5),
  etavec = c(0.25, 0.5, 0.75, 1.0),
  c0vec = c(0.25, 0.5, 0.75, 1.0),
  rvec = c(0.25, 0.5, 0.75, 1.0)
)

res <- list()
start <- 1
if (retry && file.exists(fn_raw)) {
  res <- readRDS(fn_raw)
  start <- length(res) + 1
}

if (start <= nrow(dd)) {
  
  cat("n R0 kappa eta c0 r total_pops infected_patches total_inf total_pops_qe infected_patches_qe total_inf_qe\n")
  for (i in start:nrow(dd)) {
    params <- params0
    params[["R0"]] <- dd$R0vec[i]
    params[["kappa"]] <- dd$kappavec[i]
    params[["eta"]] <- dd$etavec[i]
    params[["c0"]] <- dd$c0vec[i]
    params[["r"]] <- dd$rvec[i]
    
    res[[i]] <- mult_sim_mp(nsim = nsim, ncores = ncores, 
                            params = params, seed = 100 + i)
    cat(i, unlist(dd[i,]), sumfun1(res[[i]]), "\n")
    saveRDS(res, fn_raw)
  }
}

## summarize
res <- lapply(res, sumfun1)
dd <- data.frame(dd, do.call(rbind, res))

## attach metadata as attributes
attr(dd, "ncores") <- ncores
attr(dd, "nsim") <- nsim
attr(dd, "params0") <- params0

saveRDS(dd, fn_sum)