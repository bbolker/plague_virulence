source("simulation_funs.R")

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/raw", showWarnings = FALSE, recursive = TRUE)

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

f_ith <- function(i) file.path("outputs/raw", sprintf("res_%05d.rds", i))

cat("n R0 kappa eta c0 r mean_extinct_time extinction_rate n_extinct n_persist total_pops infected_patches total_inf\n")

for (i in 1:nrow(dd)) {
  
  if (retry && file.exists(f_ith(i))) next
  
  params <- params0
  params[["R0"]] <- dd$R0vec[i]
  params[["kappa"]] <- dd$kappavec[i]
  params[["eta"]] <- dd$etavec[i]
  params[["c0"]] <- dd$c0vec[i]
  params[["r"]] <- dd$rvec[i]
  
  one <- mult_sim_mp(nsim = nsim, ncores = ncores,
                     params = params, seed = 100 + i)
  
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