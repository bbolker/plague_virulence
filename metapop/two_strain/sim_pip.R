#!/usr/bin/env Rscript
library(optparse)
op.parser <- OptionParser(prog="sim_twostrain_poisson",
                          option_list = list(
                            make_option(c("-o", "--output"), "store", help = "output file name", default = "outputs_pip/sim_pip_alpha_rho.rds"),
                            make_option(c("-f", "--fixedpop"), "store_true", help="fix population size", default = FALSE),
                            make_option(c("-c", "--ncores"), "store", help="number of cores", default = 4L),
                            make_option(c("-n", "--nsim"), "store", help="number of simulations", default = 200L),
                            make_option(c("-r", "--rerun_existing"), "store_true", help="rerun sims for which output already exists", default = FALSE),
                            make_option("--use_r", "store_true", help="use pure-R simulator instead of C++ (slower)", default = FALSE)
                          )
                          )
opt <- parse_args(op.parser)

source("simulation_funs.R")

# Create directories for PIP results
dir.create("outputs_pip", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs_pip/raw", showWarnings = FALSE, recursive = TRUE)

fn_sum <- opt$output

ncores <- opt$ncores
nsim <- opt$nsim

# Define full parameter grid for PIPs across different alpha and rho values
dd <- expand.grid(
  R01vec = seq(1.1, 3.0, by = 0.1)
, R02vec = seq(1.1, 3.0, by = 0.1)
  ## next two are intentionally out of order, to put the 'most interesting' alpha/rho values first
, alphavec = c(5e-5, 1e-5, 1e-4)
, rhovec = c(6, 2, 8)
, c0vec = c(0.2)
, rvec = c(0.5)
)
dd$run <- seq(nrow(dd)) ## in case we re-order dd later

attr(dd, "ncores") <- ncores
attr(dd, "nsim") <- nsim
attr(dd, "params0") <- params0_2strain

saveRDS(dd, "sim_pip_design.rds")
f_ith <- function(i) file.path("outputs_pip/raw", sprintf("res_%05d.rds", i))

cat(sprintf("Starting Multi-PIP sweep with %d combinations...\n", nrow(dd)))
flush.console()

for (i in 1:nrow(dd)) {

  if (!opt$rerun_existing && file.exists(f_ith(i))) next
  
  params <- params0_2strain 
  
  params[["R01"]] <- dd$R01vec[i]
  params[["R02"]] <- dd$R02vec[i]
  params[["alpha"]] <- dd$alphavec[i]
  params[["rho"]] <- dd$rhovec[i]
  params[["c0"]] <- dd$c0vec[i]
  params[["r"]] <- dd$rvec[i]
  
  one <- mult_sim_2strain(nsim = nsim, ncores = ncores,
                          params = params, seed = 100 + i,
                          use_cpp = !opt$use_r)
  
  # Atomic write to prevent file corruption
  tmp <- paste0(f_ith(i), ".tmp")
  saveRDS(one, tmp)
  file.rename(tmp, f_ith(i))
  
  # Console progress update
  cat(sprintf("Done %d / %d (alpha=%e, rho=%.1f | R01=%.1f, R02=%.1f)\n", 
              i, nrow(dd), params[["alpha"]], params[["rho"]], params[["R01"]], params[["R02"]]))
  flush.console()
}

# Summarize results
res_list <- vector("list", nrow(dd))
for (i in 1:nrow(dd)) {
  if (file.exists(f_ith(i))) {
    res_list[[i]] <- readRDS(f_ith(i))
  }
}

res_sum <- lapply(res_list, function(x) {
  if (is.null(x)) rep(NA_real_, 10) else sumfun_2strain(x)
})

dd_out <- data.frame(dd, do.call(rbind, res_sum))

attr(dd_out, "ncores") <- ncores
attr(dd_out, "nsim") <- nsim
attr(dd_out, "params0") <- params0_2strain

saveRDS(dd_out, fn_sum)
cat("Multi-PIP sweep complete. Saved to", fn_sum, "\n")
