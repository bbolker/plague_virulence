## FIXME: add optparse stuff to make inputs/outputs flexible

source("simulation_funs.R")
fn_sum <- "outputs_pip/sim_pip_alpha_rho_tmp.rds"
dd <- readRDS("sim_pip_design.rds")
f_ith <- function(i) file.path("outputs_pip/raw", sprintf("res_%05d.rds", i))
## Summarize results
res_list <- vector("list", nrow(dd))
for (i in 1:nrow(dd)) {
  if (file.exists(f_ith(i))) {
    res_list[[i]] <- readRDS(f_ith(i))
  }
}

res_sum <- lapply(res_list, function(x) {
  if (is.null(x)) rep(NA_real_, 10) else sumfun_2strain(x)
})

res0 <- do.call(rbind, res_sum)
dd_out <- data.frame(dd, res0) |>
  dplyr::arrange(rhovec, alphavec, R01vec, R02vec)

for (p in c("ncores", "nsim", "params0")) {
  attr(dd_out, p) <- attr(dd, p)
}
saveRDS(dd_out, fn_sum)
