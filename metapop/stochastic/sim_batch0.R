source("simulation_funs.R")
library("dplyr")
library("purrr")

R0vec <- c(1.2,1.5,2.0,2.5,3)
res <- list()
for (i in seq_along(R0vec)) {
  params <- params0
  params[["R0"]] <- R0vec[i]
  res[[i]] <- mult_sim_mp(nsim = 20, ncores = 10, params = params, seed = 101+i)
  saveRDS(res, "sim_batch0_raw.rds")
}


plotfun1(res[[1]])
plotfun1(res[[3]])
plotfun1(res[[5]])
res2 <- (res
  |> map(sumfun1)
  |> setNames(R0vec)
  |> bind_rows(.id = "R0")
  |> mutate(across(R0, ~as.numeric(as.character(.))))
)

saveRDS(res2, "sim_batch0.rds")

