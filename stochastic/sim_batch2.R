#R0,alpha,K

source("simulation_funs.R")

fn_raw <- "sim_batch3_raw.rds"
fn_sum <- "sim_batch3.rds"

retry <- TRUE    ## checkpoint/pick up from previous runs?
ncores <- 16
nsim <- 120


R0vec <- seq(1.2, 3, by = 0.1)
alphavec <- 10^(-4:0)
Kvec <- 10^(2:10)

dd <- expand.grid(R0vec = R0vec,
                  alphavec = alphavec,
                  Kvec = Kvec)

res <- list()
start <- 1
if (retry && file.exists(fn_raw)) {
  res <- readRDS(fn_raw)
  start <- length(res) + 1
}

if (start <= nrow(dd)) {
  
  cat("i R0 alpha K total_pops inf_patch total_inf\n")
  for (i in start:nrow(dd)) {
    params <- params0  ## default parameters, from simulation_funs.R
    params[["R0"]] <- dd$R0vec[i]
    params[["alpha"]] <- dd$alphavec[i]
    params[["K"]] <- dd$Kvec[i]
    
    res[[i]] <- mult_sim_mp(nsim = nsim, ncores = ncores, params = params, seed = 100 + i)
    cat(i, unlist(dd[i, ]), sumfun1(res[[i]]), "\n")
    saveRDS(res, fn_raw)
  }
}

## summarize 
res_sum <- lapply(res, sumfun1)
dd_out <- data.frame(dd, do.call(rbind, res_sum))

## attach metadata as attributes
attr(dd_out, "ncores") <- ncores
attr(dd_out, "nsim") <- nsim
attr(dd_out, "params0") <- params0

saveRDS(dd_out, fn_sum)
cat("Saved summary to", fn_sum, "\n")