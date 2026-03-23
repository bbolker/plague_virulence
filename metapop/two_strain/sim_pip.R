source("sim_fun.R")

# Create directories for PIP results
dir.create("outputs_pip", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs_pip/raw", showWarnings = FALSE, recursive = TRUE)

fn_sum <- "outputs_pip/sim_pip_alpha_rho.rds"

retry <- TRUE
ncores <- 16
nsim <- 200

# Define full parameter grid for PIPs across different alpha and rho values
dd <- expand.grid(
  R01vec = seq(1.1, 3.0, by = 0.1),
  R02vec = seq(1.1, 3.0, by = 0.1),
  alphavec = c(1e-5, 5e-5, 1e-4),
  rhovec = c(2, 4, 6),
  c0vec = c(0.2),
  rvec = c(0.5)
)

f_ith <- function(i) file.path("outputs_pip/raw", sprintf("res_%05d.rds", i))

cat(sprintf("Starting Multi-PIP sweep with %d combinations...\n", nrow(dd)))
flush.console()

for (i in 1:nrow(dd)) {
  
  if (retry && file.exists(f_ith(i))) next
  
  params <- params0_2strain 
  
  params[["R01"]] <- dd$R01vec[i]
  params[["R02"]] <- dd$R02vec[i]
  params[["alpha"]] <- dd$alphavec[i]
  params[["rho"]] <- dd$rhovec[i]
  params[["c0"]] <- dd$c0vec[i]
  params[["r"]] <- dd$rvec[i]
  
  one <- mult_sim_2strain(nsim = nsim, ncores = ncores,
                          params = params, seed = 100 + i)
  
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