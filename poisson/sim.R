#!/usr/bin/env Rscript
library(optparse)
op.parser <- OptionParser(prog="sim_onestrain_poisson",
                          option_list = list(
                            make_option(c("-o", "--output"), "store", help = "output file name", default = "outputs/sim.rds"),
                            make_option(c("-f", "--fixedpop"), "store_true", help="fix population size", default = FALSE),
                            make_option(c("-c", "--ncores"), "store", help="number of cores", default = 4L),
                            make_option(c("-n", "--nsim"), "store", help="number of simulations", default = 200L),
                            make_option(c("-r", "--rerun_existing"), "store_true", help="rerun sims for which output already exists", default = FALSE)
                          )
                          )
opt <- parse_args(op.parser)

print(opt)

source("simulation_funs.R")

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/raw", showWarnings = FALSE, recursive = TRUE)

fn_sum <- opt$output

ncores <- opt$ncores
nsim <- opt$nsim
fixedpop <- opt$fixedpop

print(opt)

dd <- expand.grid(
  R0vec = seq(1.1, 3.0, by = 0.1),
  ## was: -6 to -3, but min and max are boring
  alphavec = 10^seq(-5.5, -3.5, by = 0.5),  
  rhovec = seq(2 , 10 , by=2),       
  c0vec = c(0.2),
  rvec = c(0.5)
)

f_ith <- function(i) file.path("outputs/raw", sprintf("res_%05d.rds", i))

cat("n R0 alpha rho c0 r mean_extinct_time extinction_rate n_extinct n_persist total_pops infected_patches total_inf\n")

for (i in 1:nrow(dd)) {
  
  if (!opt$rerun_existing && file.exists(f_ith(i))) next
  
  params <- params0
  params[["R0"]] <- dd$R0vec[i]
  params[["alpha"]] <- dd$alphavec[i]
  params[["rho"]] <- dd$rhovec[i]
  params[["c0"]] <- dd$c0vec[i]
  params[["r"]] <- dd$rvec[i]

  ## cat("using ", ncores, "cores\n")
  
  one <- mult_sim_mp(nsim = nsim, ncores = ncores,
                     params = params, seed = 100 + i, fixed_pop = fixedpop)
  
  # Atomic write: tmp + rename
  tmp <- paste0(f_ith(i), ".tmp")
  saveRDS(one, tmp)
  file.rename(tmp, f_ith(i))
  
  cat(i, unlist(dd[i, ]), sumfun(one), "\n")
}

# Summarize: read all existing results
res_list <- vector("list", nrow(dd))
for (i in 1:nrow(dd)) {
  if (file.exists(f_ith(i))) {
    res_list[[i]] <- readRDS(f_ith(i))
  }
}

res_sum <- lapply(res_list, function(x) {
  if (is.null(x)) rep(NA_real_, 7) else sumfun(x)
})

dd <- data.frame(dd, do.call(rbind, res_sum))

## attach metadata as attributes
attr(dd, "ncores") <- ncores
attr(dd, "nsim") <- nsim
attr(dd, "params0") <- params0

saveRDS(dd, fn_sum)
