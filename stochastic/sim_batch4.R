source("simulation_funs.R")

fn_raw <- "outputs/sim_batch4_raw.rds"
fn_sum <- "outputs/sim_batch4.rds"

retry <- TRUE  ## checkpoint/pick up from previous runs?
ncores <- 16
nsim <- 120
dd <- expand.grid(R0vec = seq(1.2, 3, by = 0.1),
                  c0vec = seq(0,1,by=0.1))


res <- list()
start <- 1
if (retry && file.exists(fn_raw)) {
  res <- readRDS(fn_raw)
  start <- length(res) + 1
}

if (start <= nrow(dd)) {
  
  cat("n R0 c0 total_pops inf_patch total_inf\n")
  for (i in start:nrow(dd)) {
    params <- params0  ## default parameters, from simulation_funs.R
    params[["R0"]] <- dd$R0vec[i]
    params[["c0"]] <- dd$c0vec[i]
    ## cat(nsim, ncores, "\n")
    ## print(params)
    res[[i]] <- mult_sim_mp(nsim = nsim, ncores = ncores, params = params, seed = 100 + i)
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

