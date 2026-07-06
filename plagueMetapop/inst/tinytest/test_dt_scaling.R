library(plagueMetapop)

## discrete_run() rescales the step column by dt when dt != 1,
## converting integer step indices to disease-generation units.

nt <- 5L

r_scaled <- discrete_run(
  n_patch   = 1L,
  nt        = nt,
  dt        = 0.2,
  def_file  = "euler_odin_def.R",
  stop_cond = NULL,
  seed      = 1L
)

r_unit <- discrete_run(
  n_patch   = 1L,
  nt        = nt,
  dt        = 1,
  def_file  = "euler_odin_def.R",
  stop_cond = NULL,
  seed      = 1L
)

expect_equal(sort(unique(r_scaled$step)), seq(0, nt * 0.2, by = 0.2),
             info = "dt=0.2: step values are scaled to disease-generation units")

expect_equal(sort(unique(r_unit$step)), seq(0L, nt),
             info = "dt=1: step values are unchanged")

expect_equal(nrow(r_scaled), nrow(r_unit),
             info = "dt scaling does not change the number of rows")
