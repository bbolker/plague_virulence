## pick some points of interest within the PIPs and see what the dynamics actually look like.
## e.g.

dd <- readRDS("sim_pip_design.rds")
par1 <- list(alphavec=5e-5, rhovec=6, R01vec=2.2, R02vec = 2.4)

df_match <- function(dd, pars, tol = 1e-8) {
  x <- lapply(names(pars),
         \(n) {abs(dd[,n] - pars[[n]])< tol}) |>
    as.data.frame() |>
    as.matrix() |>
    apply(1, all) |>
    which()
  return(x)
}

## everyone goes extinct right very quickly ... why?
i <- df_match(dd, par1)
f_ith <- function(i) file.path("outputs_pip/raw", sprintf("res_%05d.rds", i))

xx <- readRDS(f_ith(i))
str(xx$patches1)
matplot(xx$patches1, type = "l")
matplot(xx$patches2, type = "l")

## everyone is in the 'dead zone'.  Do we need to be on the shoulders???
