source("sim_fun.R")

dir.create("outputs_2strain", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs_2strain/raw", showWarnings = FALSE, recursive = TRUE)

fn_sum <- "outputs_2strain/sim_2strain.rds"

retry <- TRUE
ncores <- 16
nsim <- 200

# Same grid as single strain
dd <- expand.grid(
  R0vec = seq(1.3, 3.0, by = 0.1),
  alphavec = 10^seq(-6, -3, by = 0.5),  
  rhovec = seq(2 , 10 , by=2),       
  c0vec = c(0.2),
  rvec = c(0.5)
)

f_ith <- function(i) file.path("outputs_2strain/raw", sprintf("res_%05d.rds", i))

cat("Starting 2-strain parallel sweep...\n")

for (i in 1:nrow(dd)) {
  
  if (retry && file.exists(f_ith(i))) next
  
  params <- params0_2strain # Make sure your base 2-strain params are defined
  
  # Set R0 for Resident, and R0 - 0.2 for Invader
  params[["R01"]] <- dd$R0vec[i]
  params[["R02"]] <- dd$R0vec[i] - 0.2
  
  params[["alpha"]] <- dd$alphavec[i]
  params[["rho"]] <- dd$rhovec[i]
  params[["c0"]] <- dd$c0vec[i]
  params[["r"]] <- dd$rvec[i]
  
  # Call the 2-strain parallel wrapper
  one <- mult_sim_2strain(nsim = nsim, ncores = ncores,
                          params = params, seed = 100 + i)
  
  # Atomic write
  tmp <- paste0(f_ith(i), ".tmp")
  saveRDS(one, tmp)
  file.rename(tmp, f_ith(i))
  
  cat(sprintf("Done %d / %d (R01=%.1f, R02=%.2f)\n", i, nrow(dd), params[["R01"]], params[["R02"]]))
}

# Summarize
res_list <- vector("list", nrow(dd))
for (i in 1:nrow(dd)) {
  if (file.exists(f_ith(i))) {
    res_list[[i]] <- readRDS(f_ith(i))
  }
}

res_sum <- lapply(res_list, function(x) {
  if (is.null(x)) rep(NA_real_, 8) else sumfun_2strain(x)
})

dd_out <- data.frame(dd, do.call(rbind, res_sum))

attr(dd_out, "ncores") <- ncores
attr(dd_out, "nsim") <- nsim
attr(dd_out, "params0") <- params0_2strain

saveRDS(dd_out, fn_sum)
cat("Sweep complete. Saved to", fn_sum, "\n")